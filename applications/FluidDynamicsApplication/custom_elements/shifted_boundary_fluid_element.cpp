// Project includes
#include "includes/kratos_flags.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
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
    KRATOS_TRY;

    // Call the base element initialize method to set the constitutive law
    TBaseElement::Initialize(rCurrentProcessInfo);

    // Initialize the ELEMENTAL_DISTANCES variable (make it threadsafe)
    if (!this->Has(ELEMENTAL_DISTANCES)) {
        Vector zero_vector(NumNodes, 0.0);
        this->SetValue(ELEMENTAL_DISTANCES, zero_vector);
    }

    KRATOS_CATCH("");
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
    // Note that the INTERFACE flag is assumed to be set in the 1st layer of elements attached to the surrogate boundary e.g. by the ShiftedBoundaryMeshlessInterfaceUtility.
    // At the faces of these INTERFACE elements which are attached to BOUNDARY elements (surrogate interface gamma_tilde), the boundary flux contributions is added as a surrogate.
    if (this->Is(INTERFACE)) {
        // Initialize the element data
        ShiftedBoundaryElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        this->InitializeGeometryData(data);

        // Get the surrogate faces local IDs.
        // Note that it might happen that an INTERFACE element has no surrogate face (i.e. a unique node in the surrogate skin)
        const auto sur_bd_ids_vect = GetSurrogateFacesIds();
        if (sur_bd_ids_vect.size() != 0) {
            // Get the parent geometry data
            double size_parent;
            const auto& r_parent_geom = this->GetGeometry();
            array_1d<double, NumNodes> N_parent;
            BoundedMatrix<double, NumNodes, Dim> DN_DX_parent;  // each node's shape function derivatives at the respective node
            GeometryUtils::CalculateGeometryData(r_parent_geom, DN_DX_parent, N_parent, size_parent);
            const auto& r_boundaries = r_parent_geom.GenerateBoundariesEntities();
            DenseMatrix<unsigned int> nodes_in_faces;
            r_parent_geom.NodesInFaces(nodes_in_faces);

            // Initialize counter for the surrogate integration point
            std::size_t surrogate_pt_index = data.PositiveSideWeights.size();
            //TODO: choose - Integration method for a surrogate boundary face
            const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_1;  // r_boundaries[0].GetDefaultIntegrationMethod();

            // Loop the surrogate faces of the element
            // NOTE that there is the chance that the surrogate face is not unique
            for (std::size_t sur_bd_id : sur_bd_ids_vect) {
                // Get the current surrogate face geometry information
                const auto& r_bd_geom = r_boundaries[sur_bd_id];
                //TODO: ERROR in laplacian/small_displ?? face is column not row, for triangle row works as well because of symmetry
                const DenseVector<std::size_t> bd_local_ids = column(nodes_in_faces, sur_bd_id);  // elemental/ local node IDs of face's nodes - first entry is local ID of the face's opposite point (for tri and tetra)
                const auto& r_integration_points = r_bd_geom.IntegrationPoints(integration_method);

                // Get the gradient of the node opposite of the surrogate face
                // NOTE this gradient is used to calculate the normal as n = - DN_DX_opposite_node / norm_2(DN_DX_opposite_node)
                // NOTE this works because DN_DX of a node is calculated as the normal of the opposite face (cross product for Tetrahedra3D4N)
                BoundedVector<double,Dim> DN_DX_opposite_node = row(DN_DX_parent, bd_local_ids[0]);
                const double h_sur_bd = 1.0 / norm_2(DN_DX_opposite_node);
                BoundedVector<double,Dim> normal_sur_bd = - DN_DX_opposite_node * h_sur_bd;

                // Calculate a scaling to eliminate the change of area from global to local coordinate system of boundary's face (taken care of by element's detJ)
                // double sum_int_pt_weights = 0.0;
                // for (auto& r_int_pt : r_integration_points) {
                //    sum_int_pt_weights += r_int_pt.Weight();
                // }

                //TODO: choose weight calculation; based on element's detJ or boundary's detJ?
                // NOTE that the integration weight is calculated as detJ of the parent multiplied by the global size of the surrogate boundary face: Dim * Parent domain volume * norm(DN_DX_opposite_node)
                // NOTE that detJ is only constant inside the element for simplex elements and only for Triangle2D3N and Tetrahedra3D4N Dim*size_parent results in the correct multiplier for the norm
                // Triangle2D3N: 2*1/2*detJ*length_of_line | Tetrahedra3D4N: 3*1/6*detJ*area_of_parallelogram
                const double int_pt_weight = Dim * size_parent / h_sur_bd;
                //const double weight = Dim * size_parent / h_sur_bd / sum_int_pt_weights;

                // Get detJ for all integration points of the surrogate boundary face
                //VectorType int_pt_detJs;
                //r_bd_geom.DeterminantOfJacobian(int_pt_detJs, integration_method);

                // Loop over the integration points of the surrogate boundary face for the numerical integration of the surrogate boundary flux
                //TODO: There is only one integration point !?
                for (std::size_t i_int_pt = 0; i_int_pt < r_integration_points.size(); ++i_int_pt) {
                    // Scale integration point weight (necessary if it's more than one integration point per boundary face!)
                    // const double int_pt_weight = weight * r_integration_points[i_int_pt].Weight();
                    // Calculate integration point weight by multiplying its detJ with its weight
                    // const double int_pt_weight = int_pt_detJs[i_int_pt] * r_integration_points[i_int_pt].Weight();

                    // Compute the local coordinates of the integration point in the parent element's geometry
                    Geometry<Node>::CoordinatesArrayType aux_global_coords = ZeroVector(3);
                    r_bd_geom.GlobalCoordinates(aux_global_coords, r_integration_points[i_int_pt].Coordinates());
                    Geometry<Node>::CoordinatesArrayType int_pt_local_coords_parent = ZeroVector(3);
                    r_parent_geom.PointLocalCoordinates(int_pt_local_coords_parent, aux_global_coords);

                    // Get N of the element at the integration point
                    MatrixType int_pt_N_parent = ZeroMatrix(1, NumNodes);
                    VectorType aux_N_parent;
                    r_parent_geom.ShapeFunctionsValues(aux_N_parent, int_pt_local_coords_parent);
                    for (std::size_t i_node = 0; i_node < NumNodes; i_node++) {
                        int_pt_N_parent(0, i_node) = aux_N_parent(i_node);
                    }

                    // Get DN_DX of the element at the integration point
                    MatrixType int_pt_DN_DX_parent = ZeroMatrix(NumNodes, Dim), aux_DN_DX_parent = ZeroMatrix(NumNodes, Dim), aux_DN_DXi_parent, aux_J_parent, aux_J_inv_parent;
                    double aux_detJ_parent;
                    r_parent_geom.ShapeFunctionsLocalGradients(aux_DN_DXi_parent, int_pt_local_coords_parent);
                    r_parent_geom.Jacobian(aux_J_parent, int_pt_local_coords_parent);
                    MathUtils<double>::InvertMatrix(aux_J_parent, aux_J_inv_parent, aux_detJ_parent);
                    aux_DN_DX_parent = prod(aux_DN_DXi_parent, aux_J_inv_parent);
                    for (std::size_t d = 0; d < Dim; ++d ) {
                        for (std::size_t i_node = 0; i_node < NumNodes; ++i_node) {
                            int_pt_DN_DX_parent(i_node, d) = aux_DN_DX_parent(i_node, d);
                        }
                    }

                    // Update the element's data
                    this->UpdateIntegrationPointData(data, surrogate_pt_index++, int_pt_weight, row(int_pt_N_parent,0), int_pt_DN_DX_parent);

                    // Calculate the surrogate boundary traction as t_i = tau_ij - p_h * I_ij) * n_j with tau_ij = C:delta^s u_h taking N, DN_DX and Weight from the element's data
                    this->AddBoundaryTraction(data, normal_sur_bd, rLeftHandSideMatrix, rRightHandSideVector);
                }
            }
        }

        // ALTERNATIVE (same result) - TODO: delete?
        /*if (sur_bd_ids_vect.size() != 0) {
            // Initialize the element data
            ShiftedBoundaryElementData data;
            data.Initialize(*this, rCurrentProcessInfo);
            //this->CalculateGeometryData(data.PositiveSideWeights, data.PositiveSideN, data.PositiveSideDNDX);  //TODO maybe necessary for material response

            // Get the parent geometry data and calculate its material response (for C matrix)
            double size_parent;
            const auto& r_geom = this->GetGeometry();
            array_1d<double, NumNodes> N_parent;
            BoundedMatrix<double, NumNodes, Dim> DN_DX_parent;
            GeometryUtils::CalculateGeometryData(r_geom, DN_DX_parent, N_parent, size_parent);
            this->CalculateMaterialResponse(data);

            // Auxilary LHS matrix for summation over the element's surrogate boundaries
            BoundedMatrix<double, LocalSize, LocalSize> aux_LHS = ZeroMatrix(LocalSize, LocalSize);

            // Set strain matrix
            // Thereby, we calculate the stress at the element midpoint (not at integration point)
            // Note that in here we are assuming constant strain kinematics (which is correct for linear shape functions)
            BoundedMatrix<double, StrainSize, LocalSize> B_matrix;
            FluidElementUtilities<NumNodes>::GetStrainMatrix(DN_DX_parent, B_matrix);

            const auto &r_boundaries = r_geom.GenerateBoundariesEntities();
            DenseMatrix<unsigned int> nodes_in_faces;
            r_geom.NodesInFaces(nodes_in_faces);

            // Loop the surrogate faces
            // Note that there is the chance that the surrogate face is not unique
            for (std::size_t sur_bd_id : sur_bd_ids_vect) {
                // Get the current surrogate face geometry information
                const auto& r_sur_bd_geom = r_boundaries[sur_bd_id];
                const std::size_t n_bd_points = r_sur_bd_geom.PointsNumber();  // number of nodes of the surrogate face
                const DenseVector<std::size_t> sur_bd_local_ids = row(nodes_in_faces, sur_bd_id);
                const auto& r_sur_bd_N = r_sur_bd_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_1);  // matrix of SF values F_{ij}, where i is integration pt index and j is SF index

                // Get the gradient of the node opposite of the surrogate face
                // Note that this is used to calculate the normal as n = - DN_DX_opposite_node / norm_2(DN_DX_opposite_node)
                const BoundedVector<double,Dim> DN_DX_opposite_node = row(DN_DX_parent, sur_bd_local_ids[0]);
                const double h_sur_bd = 1.0 / norm_2(DN_DX_opposite_node);
                BoundedVector<double,Dim> normal_sur_bd = - DN_DX_opposite_node * h_sur_bd;

                // Note that the integration weight is calculated as Dim * Parent domain size * norm(DN_DX_opposite_node)
                const double weight = Dim * size_parent / h_sur_bd;

                // Compute the required projections using the boundary's normal, the elements constitutive matrix (C) and the strain matrix (B)
                BoundedMatrix<double, Dim, StrainSize> voigt_normal_projection_matrix = ZeroMatrix(Dim, StrainSize);
                FluidElementUtilities<NumNodes>::VoigtTransformForProduct(normal_sur_bd, voigt_normal_projection_matrix);
                const BoundedMatrix<double, Dim, StrainSize> aux_matrix_AC = prod(voigt_normal_projection_matrix, data.C);
                const BoundedMatrix<double, Dim, LocalSize> aux_matrix_ACB = prod(aux_matrix_AC, B_matrix);

                // Fill the shape functions auxiliary transpose matrix and the pressure to Voigt notation operator matrix for the surrogate boundary integration point
                // Note that the local face ids. are already taken into account in the assembly ??
                BoundedMatrix<double, LocalSize, Dim> N_aux_trans = ZeroMatrix(LocalSize, Dim);
                BoundedMatrix<double, StrainSize, LocalSize> pres_to_voigt_matrix_op = ZeroMatrix(StrainSize, LocalSize);
                std::size_t i_local_id;
                for (std::size_t i_bd_node = 0; i_bd_node < n_bd_points; ++i_bd_node) {
                    i_local_id = sur_bd_local_ids[i_bd_node+1];
                    for (std::size_t d = 0; d < Dim; ++d) {
                        N_aux_trans(i_local_id*BlockSize+d, d) = r_sur_bd_N(0, i_bd_node);  // TODO only one integration point ?!?
                        pres_to_voigt_matrix_op(d, i_local_id*BlockSize+Dim) = r_sur_bd_N(0, i_bd_node);
                    }
                }

                // Contribution coming from the shear stress operator
                aux_LHS += weight * prod(N_aux_trans, aux_matrix_ACB);

                // Contribution coming from the pressure terms
                const BoundedMatrix<double, LocalSize, StrainSize> N_voigt_proj_matrix = prod(N_aux_trans, voigt_normal_projection_matrix);
                aux_LHS -= weight * prod(N_voigt_proj_matrix, pres_to_voigt_matrix_op);
            }

            // Add boundary traction of the element's surrogate boundaries to the system
            array_1d<double,LocalSize> values;
            this->GetCurrentValuesVector(data,values);
            rLeftHandSideMatrix -= aux_LHS;
            rRightHandSideVector += prod(aux_LHS, values);
        }*/
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
    //TODO calculate cutted area?

    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Calculate(
    const Variable<array_1d<double, 3>> &rVariable,
    array_1d<double, 3> &rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    rOutput = ZeroVector(3);

    // If the element is split, integrate sigma.n over the interface
    // Note that in the Ausas formulation (discontinuous for thin-walled), both interface sides need to be integrated

    //TODO take care of drag force calculation see embedded element??

    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Calculate(
    const Variable<Vector> &rVariable,
    Vector& rOutput,
    const ProcessInfo &rCurrentProcessInfo)
{
    TBaseElement::Calculate(rVariable, rOutput, rCurrentProcessInfo);
}

template <class TBaseElement>
void ShiftedBoundaryFluidElement<TBaseElement>::Calculate(
    const Variable<Matrix> &rVariable,
    Matrix& rOutput,
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
        "required_variables"         : ["DISTANCE","VELOCITY","PRESSURE","MESH_VELOCITY","MESH_DISPLACEMENT"],
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
    rData.PositiveIndices.clear();
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
    this->CalculateGeometryData(rData.PositiveSideWeights, rData.PositiveSideN, rData.PositiveSideDNDX);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template<class TBaseElement>
std::vector<std::size_t> ShiftedBoundaryFluidElement<TBaseElement>::GetSurrogateFacesIds()
{
    const std::size_t n_faces = Dim + 1;
    auto& r_neigh_elems = this->GetValue(NEIGHBOUR_ELEMENTS);

    // Check the current element faces
    // Note that we rely on the fact that the neighbors are sorted according to the faces
    std::vector<std::size_t> surrogate_faces_ids;
    for (std::size_t i_face = 0; i_face < n_faces; ++i_face) {
        auto p_neigh_elem = r_neigh_elems(i_face).get();
        if (p_neigh_elem != nullptr && p_neigh_elem->Is(BOUNDARY)) {
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
