// Project includes
#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "input_output/logger.h"
#include "utilities/geometry_utilities.h"
#include "utilities/element_size_calculator.h"

#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"

// Application includes
#include "custom_elements/shifted_boundary_fluid_element.h"
#include "custom_elements/weakly_compressible_navier_stokes.h"
#include "custom_utilities/embedded_data.h"
#include "custom_utilities/weakly_compressible_navier_stokes_data.h"
#include <cstddef>
#include <ostream>
#include <string>
#include <sys/types.h>


namespace Kratos {

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< class TBaseElement >
ShiftedBoundaryFluidElement<TBaseElement>::ShiftedBoundaryFluidElement(IndexType NewId):
    TBaseElement(NewId)
{}

template< class TBaseElement >
ShiftedBoundaryFluidElement<TBaseElement>::ShiftedBoundaryFluidElement(IndexType NewId, const NodesArrayType& ThisNodes):
    TBaseElement(NewId,ThisNodes)
{}


template< class TBaseElement >
ShiftedBoundaryFluidElement<TBaseElement>::ShiftedBoundaryFluidElement(IndexType NewId, Geometry<NodeType>::Pointer pGeometry):
    TBaseElement(NewId,pGeometry)
{}

template< class TBaseElement >
ShiftedBoundaryFluidElement<TBaseElement>::ShiftedBoundaryFluidElement(IndexType NewId, Geometry<NodeType>::Pointer pGeometry, Properties::Pointer pProperties):
    TBaseElement(NewId,pGeometry,pProperties)
{}


template< class TBaseElement >
ShiftedBoundaryFluidElement<TBaseElement>::~ShiftedBoundaryFluidElement()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< class TBaseElement >
Element::Pointer ShiftedBoundaryFluidElement<TBaseElement>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ShiftedBoundaryFluidElement>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}


template< class TBaseElement >
Element::Pointer ShiftedBoundaryFluidElement<TBaseElement>::Create(
    IndexType NewId,
    Geometry<NodeType>::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ShiftedBoundaryFluidElement>(NewId, pGeom, pProperties);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    // Call the base element initialize method to set the constitutive law
    TBaseElement::Initialize(rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Add base fluid contribution (volume integration)
    TBaseElement::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    // Check if the element belongs to the surrogate interface.
    // Note that the SBM_INTERFACE flag is assumed to be set in the 1st layer of elements attached to the surrogate boundary (BOUNDARY elements) e.g. by the ShiftedBoundaryMeshlessInterfaceUtility.
    // One or multiple faces of these SBM_INTERFACE elements are attached to SBM_BOUNDARY elements.
    // These faces are therefore a part of the surrogate interface gamma_tilde, so the boundary flux contributions is added.
    if (this->Is(SBM_INTERFACE)) {
        // Initialize the element data
        ShiftedBoundaryElementData data;
        data.Initialize(*this, rCurrentProcessInfo);

        // Get the local IDs of the SBM_INTERFACE elements of the faces attached to SBM_BOUNDARY elements.
        const auto sur_bd_ids_vect = GetSurrogateFacesIds();

        // Set integration method for surrogate boundary face
        // NOTE that second order is required here
        const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;

        // Integrate over the boundary using the element's "AddBoundaryTraction" method
        // NOTE that this requires the location of the integration points of the boundary in global coordinates, for which then the local coordinates of the element are calculated
        /*if (sur_bd_ids_vect.size() != 0) {

            // Get the parent geometry data
            double size_parent;
            const auto& r_parent_geom = this->GetGeometry();
            array_1d<double, NumNodes> N_parent;
            BoundedMatrix<double, NumNodes, Dim> DN_DX_parent;  // each node's shape function derivatives at the respective node
            GeometryUtils::CalculateGeometryData(r_parent_geom, DN_DX_parent, N_parent, size_parent);
            const auto& r_boundaries = r_parent_geom.GenerateBoundariesEntities();
            DenseMatrix<unsigned int> nodes_in_faces;
            r_parent_geom.NodesInFaces(nodes_in_faces);

            // Initialize counter for the integration point index of "UpdateIntegrationPointData"
            VectorType volume_weights;
            MatrixType shape_functions;
            GeometryData::ShapeFunctionsGradientsType shape_derivatives;
            this->CalculateGeometryData(volume_weights, shape_functions, shape_derivatives);
            std::size_t surrogate_pt_index = volume_weights.size();

            // Loop the surrogate faces of the element
            // NOTE that the element might have multiple faces attached to the boundary
            for (std::size_t sur_bd_id : sur_bd_ids_vect) {
                // Get the current surrogate face's geometry information
                const auto& r_bd_geom = r_boundaries[sur_bd_id];
                // Get local node IDs of the face
                // NOTE that the column gives the local node IDs for a face and it's first entry is the local ID of the node opposite of the face (for tri and tetra)
                const DenseVector<std::size_t> bd_local_ids = column(nodes_in_faces, sur_bd_id);
                const auto& r_integration_points = r_bd_geom.IntegrationPoints(integration_method);

                // Get the gradient of the node opposite of the surrogate face to calculate the normal
                // NOTE this works because DN_DX of a node is calculated as the normal of the opposite face (cross product for Tetrahedra3D4N)
                const BoundedVector<double,Dim> DN_DX_opposite_node = row(DN_DX_parent, bd_local_ids[0]);
                const double h_sur_bd = 1.0 / norm_2(DN_DX_opposite_node);
                BoundedVector<double,3> normal_sur_bd = ZeroVector(3);
                for (std::size_t d = 0; d < Dim; ++d) {
                    normal_sur_bd(d) = - DN_DX_opposite_node(d) * h_sur_bd;
                }

                // Get detJ for all integration points of the surrogate boundary face
                VectorType int_pt_detJs;
                r_bd_geom.DeterminantOfJacobian(int_pt_detJs, integration_method);

                // Loop over the integration points of the surrogate boundary face for the numerical integration of the surrogate boundary flux
                // Therefore the local coordinates of the element are calculated from integration point position at the surrogate boundary
                // and shape function values and gradients calculated correspondingly, in order to use the elements "AddBoundaryTraction" method.
                for (std::size_t i_int_pt = 0; i_int_pt < r_integration_points.size(); ++i_int_pt) {
                    // Calculate integration point weight by multiplying its detJ with its gauss weight
                    const double int_pt_weight = int_pt_detJs[i_int_pt] * r_integration_points[i_int_pt].Weight();

                    // Compute the local coordinates of the integration point in the parent element's geometry
                    Geometry<NodeType>::CoordinatesArrayType aux_global_coords = ZeroVector(3);
                    r_bd_geom.GlobalCoordinates(aux_global_coords, r_integration_points[i_int_pt].Coordinates());
                    Geometry<NodeType>::CoordinatesArrayType int_pt_local_coords_parent = ZeroVector(3);
                    r_parent_geom.PointLocalCoordinates(int_pt_local_coords_parent, aux_global_coords);

                    // Get N of the element at the integration point
                    MatrixType int_pt_N_parent = ZeroMatrix(1, NumNodes);
                    VectorType aux_N_parent;
                    r_parent_geom.ShapeFunctionsValues(aux_N_parent, int_pt_local_coords_parent);
                    for (std::size_t i_node = 0; i_node < NumNodes; i_node++) {
                        int_pt_N_parent(0, i_node) = aux_N_parent(i_node);
                    }

                    // Get DN_DX of the element at the integration point
                    MatrixType int_pt_DN_DX_parent = ZeroMatrix(NumNodes, Dim), aux_DN_DXi_parent, aux_J_parent, aux_J_inv_parent;
                    double aux_detJ_parent;
                    r_parent_geom.ShapeFunctionsLocalGradients(aux_DN_DXi_parent, int_pt_local_coords_parent);
                    r_parent_geom.Jacobian(aux_J_parent, int_pt_local_coords_parent);
                    MathUtils<double>::InvertMatrix(aux_J_parent, aux_J_inv_parent, aux_detJ_parent);
                    int_pt_DN_DX_parent = prod(aux_DN_DXi_parent, aux_J_inv_parent);

                    // Update the element's data
                    this->UpdateIntegrationPointData(data, surrogate_pt_index++, int_pt_weight, row(int_pt_N_parent,0), int_pt_DN_DX_parent);

                    // Calculate the surrogate boundary traction as t_i = (tau_ij - p_h * I_ij) * n_j with tau_ij = C:delta^s u_h, taking N, DN_DX and Weight from the element's data
                    this->AddBoundaryTraction(data, normal_sur_bd, rLeftHandSideMatrix, rRightHandSideVector);
                }
            }
        }*/

        // Integrate over the boundary using the faces' geometry
        // NOTE that this is the safer alternative because integration points are used directly and do not need to be found inside the element
        //TODO SPEED UP
        if (sur_bd_ids_vect.size() != 0) {

            // Get the parent geometry data. It is calculated at the element midpoint.
            double size_parent;
            const auto& r_geom = this->GetGeometry();
            array_1d<double, NumNodes> N_parent;
            BoundedMatrix<double, NumNodes, Dim> DN_DX_parent;
            GeometryUtils::CalculateGeometryData(r_geom, DN_DX_parent, N_parent, size_parent);

            // Calculate material response in order to update C matrix (constant)
            this->CalculateMaterialResponse(data);

            // Auxilary LHS matrix for summation over the element's surrogate boundaries
            BoundedMatrix<double, LocalSize, LocalSize> aux_LHS = ZeroMatrix(LocalSize, LocalSize);

            // Get strain matrix
            // NOTE that here we are assuming constant strain kinematics (which is correct for the linear shape functions of simplex elements)
            BoundedMatrix<double, StrainSize, LocalSize> B_matrix;
            FluidElementUtilities<NumNodes>::GetStrainMatrix(DN_DX_parent, B_matrix);

            const auto &r_boundaries = r_geom.GenerateBoundariesEntities();
            DenseMatrix<unsigned int> nodes_in_faces;
            r_geom.NodesInFaces(nodes_in_faces);

            // Loop the surrogate faces
            // NOTE that more than one face of the element might be part of the surrogate boundary
            for (std::size_t sur_bd_id : sur_bd_ids_vect) {
                // Get the surrogate face's geometry information
                const auto& r_sur_bd_geom = r_boundaries[sur_bd_id];
                const std::size_t n_bd_points = r_sur_bd_geom.PointsNumber();  // number of nodes of the surrogate face
                const DenseVector<std::size_t> sur_bd_local_ids = column(nodes_in_faces, sur_bd_id);

                // Get integration points and their global shape function values of the surrogate boundary face
                const auto& r_sur_bd_N = r_sur_bd_geom.ShapeFunctionsValues(integration_method);  // matrix of SF values F_{ij}, where i is integration pt index and j is SF index
                const std::vector<IntegrationPoint<3>> & r_integration_points = r_sur_bd_geom.IntegrationPoints(integration_method);

                // Get values of determinant of Jacobian for all integration points of the surrogate boundary face
                VectorType sur_bd_detJ;
                r_sur_bd_geom.DeterminantOfJacobian(sur_bd_detJ, integration_method);

                // Get the gradient of the node opposite of the surrogate face
                // NOTE that the gradient is used to calculate the normal as n = - DN_DX_opposite_node / norm_2(DN_DX_opposite_node)
                const BoundedVector<double,Dim> DN_DX_opposite_node = row(DN_DX_parent, sur_bd_local_ids[0]);
                const double h_sur_bd = 1.0 / norm_2(DN_DX_opposite_node);
                BoundedVector<double,3> normal_sur_bd = ZeroVector(3);
                for (std::size_t d = 0; d < Dim; ++d) {
                    normal_sur_bd(d) = - DN_DX_opposite_node(d) * h_sur_bd;
                }

                // Get the normal of the true boundary (Gamma) //TODO
                BoundedVector<double,3> normal_Gamma = ZeroVector(3);
                BoundedVector<double,3> center_Gamma = ZeroVector(3);
                // center_Gamma(0) = 2.75;
                // center_Gamma(2) = 1.0;
                double dot_product_Gamma_bd1 = r_sur_bd_geom[sur_bd_local_ids[1]][0] * 2.75 + r_sur_bd_geom[sur_bd_local_ids[1]][2];
                if (dot_product_Gamma_bd1 > 1.0) {
                    normal_Gamma(0) = -0.044400614;
                    normal_Gamma(2) = -0.999013806;
                } else {
                    normal_Gamma(0) =  0.044400614;
                    normal_Gamma(2) =  0.999013806;
                }

                // Compute the normal projection matrix for the local surrogate boundary
                BoundedMatrix<double, Dim, Dim> normal_proj_matrix = ZeroMatrix(Dim, Dim);
                for (std::size_t d = 0; d < Dim; ++d) {
                   normal_proj_matrix(d, d) = normal_sur_bd(d);
                }

                // Get the normal projection matrix for the local surrogate boundary in Voigt notation (A_matrix)
                BoundedMatrix<double, Dim, StrainSize> voigt_normal_projection_matrix = ZeroMatrix(Dim, StrainSize);
                FluidElementUtilities<NumNodes>::VoigtTransformForProduct(normal_sur_bd, voigt_normal_projection_matrix);

                // Get the normal projection matrix for the true boundary in Voigt notation (A_matrix_Gamma)
                BoundedMatrix<double, Dim, StrainSize> voigt_normal_proj_matrix_Gamma;
                FluidElementUtilities<3>::VoigtTransformForProduct(normal_Gamma, voigt_normal_proj_matrix_Gamma);

                // Get the tangential projection matrix (I - n x n) for the true boundary
                BoundedMatrix<double, Dim, Dim> tang_proj_matrix_Gamma;
                FluidElementUtilities<3>::SetTangentialProjectionMatrix(normal_Gamma, tang_proj_matrix_Gamma);

                // Compute auxilary matrices for the shear stress using the boundary's normal, the elements constitutive matrix (C) and the strain matrix (B)
                const BoundedMatrix<double, StrainSize, LocalSize> aux_matrix_CB = prod(data.C, B_matrix);
                //const BoundedMatrix<double, Dim, StrainSize> aux_matrix_AC = prod(voigt_normal_projection_matrix, data.C);
                const BoundedMatrix<double, Dim, LocalSize> aux_matrix_ACB = prod(voigt_normal_projection_matrix, aux_matrix_CB);  //prod(aux_matrix_AC, B_matrix);
                const BoundedMatrix<double, Dim, LocalSize> aux_matrix_ACB_Gamma = prod(voigt_normal_proj_matrix_Gamma, aux_matrix_CB);

                // Loop over the integration points of the surrogate boundary face for the numerical integration of the surrogate boundary flux
                for (std::size_t i_int_pt = 0; i_int_pt < r_integration_points.size(); ++i_int_pt) {
                    // Fill the shape functions auxiliary transpose matrix and the pressure to Voigt notation operator matrix for the surrogate boundary integration point
                    // NOTE that the local face IDs are already taken into account in the assembly
                    BoundedMatrix<double, LocalSize, Dim> N_aux_u_trans = ZeroMatrix(LocalSize, Dim);
                    //BoundedMatrix<double, StrainSize, LocalSize> pres_to_voigt_matrix_op = ZeroMatrix(StrainSize, LocalSize);
                    BoundedMatrix<double, LocalSize, Dim> N_aux_p = ZeroMatrix(Dim, LocalSize);
                    std::size_t node_id_in_element;
                    for (std::size_t i_bd_node = 0; i_bd_node < n_bd_points; ++i_bd_node) {
                        node_id_in_element = sur_bd_local_ids[i_bd_node+1];
                        for (std::size_t d = 0; d < Dim; ++d) {
                            N_aux_u_trans(node_id_in_element*BlockSize+d, d)             = r_sur_bd_N(i_int_pt, i_bd_node);
                            //pres_to_voigt_matrix_op(d, node_id_in_element*BlockSize+Dim) = r_sur_bd_N(i_int_pt, i_bd_node);
                            N_aux_p(d, node_id_in_element*BlockSize+Dim)                 = r_sur_bd_N(i_int_pt, i_bd_node);
                        }
                    }

                    // Calculate integration point weight by multiplying its detJ with its gauss weight
                    const double int_pt_weight = sur_bd_detJ[i_int_pt] * r_integration_points[i_int_pt].Weight();

                    // Contribution coming from the shear stress operator
                    aux_LHS -= int_pt_weight * prod(N_aux_u_trans, aux_matrix_ACB);

                    // Contribution coming from the pressure term
                    //const BoundedMatrix<double, LocalSize, StrainSize> N_voigt_proj_matrix = prod(N_aux_u_trans, voigt_normal_projection_matrix);
                    //aux_LHS += int_pt_weight * prod(N_voigt_proj_matrix, pres_to_voigt_matrix_op);
                    // Alternative normal projection
                    const BoundedMatrix<double, LocalSize, Dim> N_normal_proj_matrix = prod(N_aux_u_trans, normal_proj_matrix);
                    aux_LHS += int_pt_weight * prod(N_normal_proj_matrix, N_aux_p);

                    // Nitsche slip tangential penalty - Contribution coming from the traction vector tangential component (no contribution for no slip)
                    const Matrix aux_matrix_N_trans_tang_Gamma = prod(N_aux_u_trans, tang_proj_matrix_Gamma);
                    const double pen_coeff_tang_shear = 1.0;
                    aux_LHS += pen_coeff_tang_shear  * int_pt_weight * prod(aux_matrix_N_trans_tang_Gamma, aux_matrix_ACB_Gamma);
                }
            }

            // Add boundary traction of the element's surrogate boundaries to the system
            array_1d<double,LocalSize> values;
            this->GetCurrentValuesVector(data, values);
            rLeftHandSideMatrix  += aux_LHS;
            rRightHandSideVector -= prod(aux_LHS, values);
        }
    }

    KRATOS_CATCH("")
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    VectorType rRightHandSideVector;
    this->CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType rLeftHandSideMatrix;
    this->CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Calculate(
    const Variable<double> &rVariable,
    double& rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Calculate(
    const Variable<array_1d<double, 3>> &rVariable,
    array_1d<double, 3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    rOutput = ZeroVector(3);
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Calculate(
    const Variable<VectorType> &rVariable,
    VectorType& rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Calculate(
    const Variable<MatrixType> &rVariable,
    MatrixType& rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <class TBaseElement>
int ShiftedBoundaryFluidElement<TBaseElement>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    int out = ShiftedBoundaryElementData::Check(*this, rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Something is wrong with the elemental data of Element "
        << this->Info() << std::endl;

    return TBaseElement::Check(rCurrentProcessInfo);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <class TBaseElement>
const Parameters ShiftedBoundaryFluidElement<TBaseElement>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["implicit"],
        "framework"                  : "ale",
        "symmetric_lhs"              : false,
        "positive_definite_lhs"      : true,
        "output"                     : {
            "gauss_point"            : [],
            "nodal_historical"       : ["VELOCITY","PRESSURE"],
            "nodal_non_historical"   : ["EMBEDDED_VELOCITY"],
            "entity"                 : []
        },
        "required_variables"         : ["VELOCITY","PRESSURE","MESH_VELOCITY","MESH_DISPLACEMENT"],  //NOTE necessary for (dis-)continuous level set: "DISTANCE"
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3","Tetrahedra3D4"],
        "element_integrates_in_time" : true,
        "compatible_constitutive_laws": {
            "type"        : ["Newtonian2DLaw","Newtonian3DLaw","NewtonianTemperatureDependent2DLaw","NewtonianTemperatureDependent3DLaw","Euler2DLaw","Euler3DLaw"],
            "dimension"   : ["2D","3D"],
            "strain_size" : [3,6]
        },
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   :
            "This element is based on the Shifted-Boundary Method using MLS shape functions. Therfore, it adds surrogate boundary flux contributions for elements marked with the INTERFACE flag. The element is meant to be used together with the NavierStokesShiftedBoundarySolver, which sets up the ShiftedBoundaryMeshlessInterfaceUtility creating ShiftedBoundaryWallConditions for the interface."
    })");

    if (Dim == 2) {
        std::vector<std::string> dofs_2d({"VELOCITY_X","VELOCITY_Y","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"VELOCITY_X","VELOCITY_Y","VELOCITY_Z","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

template <class TBaseElement>
std::string ShiftedBoundaryFluidElement<TBaseElement>::Info() const
{
    std::stringstream buffer;
    buffer << "ShiftedBoundaryFluidElement #" << this->Id();
    return buffer.str();
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "ShiftedBoundaryFluidElement" << Dim << "D" << NumNodes << "N"
             << std::endl
             << "on top of ";
    TBaseElement::PrintInfo(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Protected functions
///////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////
// Operations

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::InitializeGeometryData(ShiftedBoundaryElementData& rData) const
{
    //NOTE necessary for (dis-)continuous level set with EmbeddedData
    /*rData.PositiveIndices.clear();
    rData.NegativeIndices.clear();

    // Number of positive and negative distance function values
    for (std::size_t i = 0; i < ShiftedBoundaryElementData::NumNodes; ++i){
        if (rData.Distance[i] > 0.0) {
            rData.NumPositiveNodes++;
            rData.PositiveIndices.push_back(i);
        }
        else {
            rData.NumNegativeNodes++;
            rData.NegativeIndices.push_back(i);
        }
    }
    rData.NumNegativeNodes = 0;
    rData.NumPositiveNodes = NumNodes;
    this->CalculateGeometryData(rData.PositiveSideWeights, rData.PositiveSideN, rData.PositiveSideDNDX);*/
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template<class TBaseElement>
std::vector<std::size_t> ShiftedBoundaryFluidElement<TBaseElement>::GetSurrogateFacesIds()
{
    const std::size_t n_faces = Dim + 1;  // NOTE this is only valid for tri and tetra
    auto& r_neigh_elems = this->GetValue(NEIGHBOUR_ELEMENTS);

    // Check the current element faces
    // NOTE that we rely on the fact that the neighbors are sorted according to the order of faces in NEIGHBOUR_ELEMENTS
    std::vector<std::size_t> surrogate_faces_ids;
    for (std::size_t i_face = 0; i_face < n_faces; ++i_face) {
        auto p_neigh_elem = r_neigh_elems(i_face).get();
        if (p_neigh_elem != nullptr && p_neigh_elem->Is(SBM_BOUNDARY)) {
            surrogate_faces_ids.push_back(i_face);
        }
    }

    return surrogate_faces_ids;
}

// serializer

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TBaseElement);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TBaseElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class ShiftedBoundaryFluidElement< WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<2,3> > >;
template class ShiftedBoundaryFluidElement< WeaklyCompressibleNavierStokes< WeaklyCompressibleNavierStokesData<3,4> > >;

///////////////////////////////////////////////////////////////////////////////////////////////////

}  // namespace Kratos
