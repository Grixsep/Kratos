//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//

// System includes
// see the header file

// External includes
// see the header file

// Include Base headers
// see the header file

// Project includes
#include "elliptic_distance_reinitialization_element.h"
#include "droplet_dynamics_application_variables.h"
#include "utilities/element_size_calculator.h"
#include "custom_utilities/custom_tetrahedra_3d_4_modified_shape_functions.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * Constructor.
 */
EDReinitializationElement::EDReinitializationElement(IndexType NewId)
    : Element(NewId)
{
}

/**
 * Constructor using an array of nodes
 */
EDReinitializationElement::EDReinitializationElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
{
}

/**
 * Constructor using Geometry
 */
EDReinitializationElement::EDReinitializationElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
}

/**
 * Constructor using Properties
 */
EDReinitializationElement::EDReinitializationElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

/**
 * Copy Constructor
 */
EDReinitializationElement::EDReinitializationElement(EDReinitializationElement const& rOther)
    : Element(rOther)
{
}

/**
 * Destructor
 */
EDReinitializationElement::~EDReinitializationElement()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
EDReinitializationElement & EDReinitializationElement::operator=(EDReinitializationElement const& rOther)
{
    BaseType::operator=(rOther);
    Flags::operator =(rOther);
    // mpProperties = rOther.mpProperties;
    return *this;
}

///@}
///@name Operations
///@{

/**
 * ELEMENTS inherited from this class have to implement next
 * Create and Clone methods: MANDATORY
 */

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer EDReinitializationElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EDReinitializationElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer EDReinitializationElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EDReinitializationElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer EDReinitializationElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EDReinitializationElement>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void EDReinitializationElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
{
    const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;

    const int num_dof = num_nodes; //(num_dim + 1)*num_nodes;
    if (rResult.size() != num_dof){
        rResult.resize(num_dof, false);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        rResult[i]  =  this->GetGeometry()[i].GetDof(DISTANCE_AUX2).EquationId();
    }
}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
void EDReinitializationElement::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const
{
    const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;

    const int num_dof = num_nodes; //(num_dim + 1)*num_nodes;
    if (rElementalDofList.size() != num_dof){
        rElementalDofList.resize(num_dof);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        rElementalDofList[i] = this->GetGeometry()[i].pGetDof(DISTANCE_AUX2);
    }
}

/**
 * ELEMENTS inherited from this class have to implement next
 * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
 * they can be managed internally with a private method to do the same calculations
 * only once: MANDATORY
 */

/**
 * this is called during the assembling process in order
 * to calculate all elemental contributions to the global system
 * matrix and the right hand side
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void EDReinitializationElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const int num_dim  = 3;
    const int num_nodes  = num_dim + 1;
    const int num_dof = num_nodes;

    const unsigned int num_faces = num_nodes; //for simplex elements
    const unsigned int num_face_nodes = num_nodes - 1;

    const double element_size = ElementSizeCalculator<3,4>::AverageElementSize(this->GetGeometry()); //MinimumElementSize!!?
    const double tolerance = 1.0e-3*element_size;

    GeometryData::ShapeFunctionsGradientsType DN_DX;
    Matrix N;
    Vector DetJ;
    Vector weights;

    Vector distances0(num_nodes);
    Vector values(num_dof);
    VectorType grad_phi_avg;    //Recovered (lumped-mass) gradient
    //VectorType grad_phi;      //Based on DN_DX*Phi
    VectorType grad_phi_old;    //Based on the unmodified DISTANCE, Only for normal calculation

    BoundedMatrix<double,num_dof,num_dof> lhs = ZeroMatrix(num_dof,num_dof);
    BoundedMatrix<double,num_dof,num_dof> lhs_Newton_Raphson = ZeroMatrix(num_dof,num_dof);
    Vector rhs = ZeroVector(num_dof);

    const GeometryData::IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_2;
    GeometryType::Pointer p_geometry = this->pGetGeometry();
    GeometryType::GeometriesArrayType faces = p_geometry->Faces();
    const unsigned int number_of_gauss_points = p_geometry->IntegrationPointsNumber(integration_method);

    // Getting data for the given geometry
    p_geometry->ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, integration_method);

    if (N.size1() != number_of_gauss_points || N.size2() != num_nodes) {
        N.resize(number_of_gauss_points,num_nodes,false);
    }
    N = p_geometry->ShapeFunctionsValues(integration_method);

    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = p_geometry->IntegrationPoints(integration_method);
    if (weights.size() != number_of_gauss_points) {
        weights.resize(number_of_gauss_points,false);
    }
    for (unsigned int gp = 0; gp < number_of_gauss_points; gp++){
        weights[gp] = DetJ[gp] * IntegrationPoints[gp].Weight();
    }

    unsigned int nneg=0, npos=0;

    double mean_curvature = 0.0;
    double mean_norm_distance_grad = 0.0;

    for (unsigned int i_node=0; i_node < num_nodes; ++i_node){
        const double dist0 = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE);
        distances0(i_node) = dist0;

        mean_curvature += (*p_geometry)[i_node].GetValue(CURVATURE);//FastGetSolutionStepValue(CURVATURE);
        mean_norm_distance_grad += norm_2((*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE_GRADIENT));

        if (dist0 > 0.0) npos += 1;
        else nneg += 1;

        const double dist_dof = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE_AUX2);
        values(i_node) = dist_dof;
    }

    mean_curvature /= static_cast<double>(num_nodes);
    //mean_curvature /= mean_norm_distance_grad;

    const double scale = 1.0e0; // For very small (micrometric) elements

    const double penalty_phi0 = scale*1.0e4/element_size; // For Nitsche's method we need 1/h

    if(mean_curvature > 0.5/element_size)   // Sharp corners
        mean_curvature = 0.5/element_size;

    const double source_coeff = 1.0*scale*mean_curvature;

    if(rLeftHandSideMatrix.size1() != num_dof)
        rLeftHandSideMatrix.resize(num_dof,num_dof,false); //resizing the system in case it does not have the right size

    if(rRightHandSideVector.size() != num_dof)
        rRightHandSideVector.resize(num_dof,false);

    const unsigned int step = rCurrentProcessInfo[FRACTIONAL_STEP];

    double diffusion = 0.0;
    double diffusion_prime_to_s = 0.0;

    for (unsigned int gp = 0; gp < number_of_gauss_points; gp++){
        grad_phi_avg = ZeroVector(num_dim);
        for (unsigned int i_node = 0; i_node < num_nodes; i_node++){
            grad_phi_avg += GetGeometry()[i_node].GetValue(DISTANCE_GRADIENT)*N(gp, i_node);
        }
        const double norm_grad_phi_avg = norm_2(grad_phi_avg);

        //grad_phi = prod(trans(DN_DX[gp]),values);
        //const double norm_grad_phi = norm_2(grad_phi);

        grad_phi_old = prod(trans(DN_DX[gp]),distances0);

        if (norm_grad_phi_avg >= 1.0) { //tolerance){
            diffusion = 1.0/norm_grad_phi_avg;
            diffusion_prime_to_s = -1.0/(norm_grad_phi_avg*norm_grad_phi_avg*norm_grad_phi_avg);
        } else{
            //diffusion = 1.2/(norm_grad_phi_avg+tolerance);
            //diffusion = 2.0-norm_grad_phi_avg;
            diffusion = 2.0 - 1.0/(2.0 - norm_grad_phi_avg);
            //diffusion_prime_to_s = -1.2/(norm_grad_phi_avg*norm_grad_phi_avg*norm_grad_phi_avg+tolerance);
            //diffusion_prime_to_s = -1.0/(norm_grad_phi_avg+tolerance);
            diffusion_prime_to_s = 1.0/(norm_grad_phi_avg+tolerance)/(2.0 - norm_grad_phi_avg)/(2.0 - norm_grad_phi_avg);
        }

        /* if (norm_grad_phi_avg > 1.0){
            diffusion = 1.0/norm_grad_phi_avg;
            diffusion_prime_to_s = -1.0/(norm_grad_phi_avg*norm_grad_phi_avg*norm_grad_phi_avg);
        } else{
            diffusion = 1.0 - 2.0*norm_grad_phi_avg*norm_grad_phi_avg*(1.0 - norm_grad_phi_avg)*(1.0 - norm_grad_phi_avg)
                + (1.0 - norm_grad_phi_avg)*norm_grad_phi_avg*norm_grad_phi_avg*norm_grad_phi_avg; //(3.0 - 2.0*norm_grad_phi)*norm_grad_phi;
            diffusion_prime_to_s = -(4.0*(1.0 - norm_grad_phi_avg)*(1.0 - norm_grad_phi_avg)
                + 7.0*norm_grad_phi_avg*(norm_grad_phi_avg - 1.0) + norm_grad_phi_avg*norm_grad_phi_avg );
        } */

        for (unsigned int i_node = 0; i_node < num_nodes; i_node++){

            double grad_Ni_dot_grad_phi = 0.0;

            for (unsigned int k_dim = 0; k_dim < num_dim; k_dim++){
                grad_Ni_dot_grad_phi += (DN_DX[gp])(i_node, k_dim) * grad_phi_avg[k_dim];
            }

            for (unsigned int j_node = 0; j_node < num_nodes; j_node++){

                double grad_Nj_dot_grad_phi = 0.0;

                for (unsigned int k_dim = 0; k_dim < num_dim; k_dim++){

                    grad_Nj_dot_grad_phi += (DN_DX[gp])(j_node, k_dim) * grad_phi_avg[k_dim];

                    if (step > 1){  // Elliptic redistancing
                        // move to the LHS for Newton-Raphson strategy:
                        rhs[i_node] += diffusion * weights(gp) * (DN_DX[gp])(i_node, k_dim) * N(gp, j_node) * grad_phi_avg[k_dim];
                        lhs(i_node, j_node) += 1.0/scale * weights(gp) * (DN_DX[gp])(i_node, k_dim) * (DN_DX[gp])(j_node, k_dim);
                        //lhs(i_node, j_node) += - diffusion * weights(gp) * (DN_DX[gp])(i_node, k_dim) * (DN_DX[gp])(j_node, k_dim);
                    }
                    else{    // Elliptic reinitialization
                        lhs(i_node, j_node) += weights(gp) * (DN_DX[gp])(i_node, k_dim) * (DN_DX[gp])(j_node, k_dim);
                    }
                }
                /* if (step > 1){
                    lhs_Newton_Raphson(i_node, j_node) -= diffusion_prime_to_s * weights(gp) * grad_Ni_dot_grad_phi * grad_Nj_dot_grad_phi;
                } */
            }
            if (step <= 1){ // For cases that source_coeff is the same for both the positive and negative domains.
                rhs[i_node] -= source_coeff * weights(gp) * N(gp, i_node);
            }
        }
    }

    if (npos != 0 && nneg != 0)
    {   
        std::vector<int> structure_node_id;
        for (unsigned int i_node = 0; i_node < num_nodes; i_node++){
            if ((*p_geometry)[i_node].GetValue(IS_STRUCTURE) == 1.0 )
                structure_node_id.push_back(i_node);
        }
        // KRATOS_INFO("EllipticDistanceReinitialization") << "this is a splited element!" << std::endl;
        CustomModifiedShapeFunctions::Pointer p_modified_sh_func =
            Kratos::make_shared<CustomTetrahedra3D4ModifiedShapeFunctions>(p_geometry, distances0, structure_node_id);

       /*  if (step <= 1){  // Uncomment if different sources are needed for the positive and negatice sides
            Matrix neg_N, pos_N;
            GeometryType::ShapeFunctionsGradientsType neg_DN_DX, pos_DN_DX;
            Vector neg_weights, pos_weights;

            p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
                pos_N,        // N
                pos_DN_DX,    // DN_DX
                pos_weights,  // weight * detJ
                integration_method);

            const std::size_t number_of_pos_gauss_points = pos_weights.size();

            for (unsigned int pos_gp = 0; pos_gp < number_of_pos_gauss_points; pos_gp++){
                for (unsigned int i_node = 0; i_node < num_nodes; i_node++){
                    // Basically, there is no need for the positive source term
                    rhs(i_node) -= 1.0e0 * source_coeff * pos_weights(pos_gp) * pos_N(pos_gp, i_node);
                }
            }

            p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
                neg_N,        // N
                neg_DN_DX,    // DN_DX
                neg_weights,  // weight * detJ
                integration_method);

            const std::size_t number_of_neg_gauss_points = neg_weights.size();

            for (unsigned int neg_gp = 0; neg_gp < number_of_neg_gauss_points; neg_gp++){
                for (unsigned int i_node = 0; i_node < num_nodes; i_node++){
                    rhs(i_node) -= 1.0e0 * source_coeff * neg_weights(neg_gp) * neg_N(neg_gp, i_node);
                }
            }
        } */

        // **********************************************************
        // **********************************************************
        // Penalizing the location of the interface
        Matrix int_N;
        GeometryType::ShapeFunctionsGradientsType int_DN_DX;
        Vector int_weights;

        p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
            int_N,
            int_DN_DX,
            int_weights,
            integration_method);

        const std::size_t number_of_int_gauss_points = int_weights.size();

        const VectorType normal0 = 1.0/norm_2(grad_phi_old) * grad_phi_old; // For a linear element, DN_DX is constant.

        for (unsigned int int_gp = 0; int_gp < number_of_int_gauss_points; int_gp++){

            for (unsigned int i_node = 0; i_node < num_nodes; i_node++){

                double normal_dot_grad_int_Ni = 0.0;

                for (unsigned int k_dim = 0; k_dim < num_dim; k_dim++){
                    normal_dot_grad_int_Ni += (int_DN_DX[int_gp])(i_node, k_dim) * normal0[k_dim];
                }

                for (unsigned int j_node = 0; j_node < num_nodes; j_node++){

                    double normal_dot_grad_int_Nj = 0.0;
                    for (unsigned int k_dim = 0; k_dim < num_dim; k_dim++){
                        normal_dot_grad_int_Nj += (int_DN_DX[int_gp])(j_node, k_dim) * normal0[k_dim];
                    }

                    lhs(i_node, j_node) += penalty_phi0*int_weights(int_gp)*int_N(int_gp, i_node)*int_N(int_gp, j_node)
                        - int_weights(int_gp)*int_N(int_gp, i_node)*normal_dot_grad_int_Nj
                        - int_weights(int_gp)*normal_dot_grad_int_Ni*int_N(int_gp, j_node);
                }
            }
        }

        // **********************************************************
        // **********************************************************
        // Penalizing the location of the contact-line
        std::vector<unsigned int> contact_line_faces;
        std::vector<unsigned int> contact_line_indices;

        std::vector<MatrixType> contact_N;                                   //std::vector for multiple contact lines
        std::vector<GeometryType::ShapeFunctionsGradientsType> contact_DN_DX;//std::vector for multiple contact lines
        std::vector<Kratos::Vector> contact_weights;                         //std::vector for multiple contact lines

        auto& neighbour_elems = this->GetValue(NEIGHBOUR_ELEMENTS);

        std::vector<Kratos::Vector> ContactTangentialsNeg; //Dummy, not needed here
        p_modified_sh_func->ComputeNegativeSideContactLineVector(contact_line_faces, ContactTangentialsNeg);

        for (unsigned int i_cl = 0; i_cl < contact_line_faces.size(); i_cl++){

            /* if (neighbour_elems[ contact_line_faces[i_cl] ].Id() == this->Id() ){ */

            //GeometryType& r_face = faces[ contact_line_faces[i_cl] ];

            unsigned int num_structure_nodes_on_contact_face = 0;

            for (unsigned int i_node = 0; i_node < num_nodes; i_node++){
                if ( (*p_geometry)[i_node].Is(BOUNDARY) && i_node != contact_line_faces[i_cl] ){//r_face[i_face_node].GetValue(IS_STRUCTURE) == 1.0 ){
                    num_structure_nodes_on_contact_face++;
                }
            }

            if (num_structure_nodes_on_contact_face == num_face_nodes){
                contact_line_indices.push_back(i_cl);}
        }
        if (contact_line_indices.size() > 0){
            KRATOS_INFO("EllipticDistanceReinitialization") << "this is a splited element!" << contact_line_indices.size() << std::endl;}

        // Call the Contact Line negative side shape functions calculator
        p_modified_sh_func->ComputeContactLineNegativeSideShapeFunctionsAndGradientsValues(
            contact_line_indices, //ADDED
            contact_N,
            contact_DN_DX,
            contact_weights,
            GeometryData::IntegrationMethod::GI_GAUSS_2);

        for (unsigned int i_cl = 0; i_cl < contact_weights.size(); i_cl++){

            const std::size_t number_of_contact_gauss_points = (contact_weights[i_cl]).size();
            for (unsigned int contact_gp = 0; contact_gp < number_of_contact_gauss_points; contact_gp++){

                for (unsigned int i_node = 0; i_node < num_nodes; i_node++){

                    double normal_dot_grad_contact_Ni = 0.0;

                    for (unsigned int k_dim = 0; k_dim < num_dim; k_dim++){
                        normal_dot_grad_contact_Ni += ( (contact_DN_DX[i_cl] )[contact_gp])(i_node, k_dim) * normal0[k_dim];
                    }

                    for (unsigned int j_node = 0; j_node < num_nodes; j_node++){

                        double normal_dot_grad_contact_Nj = 0.0;

                        for (unsigned int k_dim = 0; k_dim < num_dim; k_dim++){
                            normal_dot_grad_contact_Nj += ( (contact_DN_DX[i_cl] )[contact_gp])(j_node, k_dim) * normal0[k_dim];
                        }

                        lhs(i_node, j_node) += 1.0e0/* *element_size */*penalty_phi0*(contact_weights[i_cl])(contact_gp)*(contact_N[i_cl])(contact_gp, i_node)*(contact_N[i_cl])(contact_gp, j_node)
                            - (contact_weights[i_cl])(contact_gp)*(contact_N[i_cl])(contact_gp, i_node)*normal_dot_grad_contact_Nj
                            - (contact_weights[i_cl])(contact_gp)*normal_dot_grad_contact_Ni*(contact_N[i_cl])(contact_gp, j_node);
                    }
                }
            }
        }

    }

    /* else if (step <= 1){    // Uncomment if different sources are needed for the positive and negatice sides
        double source;
        if (npos != 0)
            source = -1.0e0; // Basically, there is no need to add positive source term
        else
            source = -1.0e0;
        for (unsigned int gp = 0; gp < number_of_gauss_points; gp++){
            for (unsigned int i_node = 0; i_node < num_nodes; i_node++){
                rhs(i_node) += source_coeff * source * weights(gp) * N(gp, i_node);
            }
        }
    } */

    // **********************************************************
    // **********************************************************
    // Imposing the Neumann BC
    if (step <= 1){
        unsigned int i_face = 0;

        while (i_face < num_faces) {

            GeometryType& r_face = faces[i_face];
            unsigned int boundary_face_node = 0;

            double contact_angle = 0.0;
            double contact_angle_weight = 0.0;
            Vector solid_normal = ZeroVector(num_dim);

            for (unsigned int i=0; i < num_face_nodes; ++i){
                if ( r_face[i].Is(BOUNDARY) || r_face[i].GetValue(IS_STRUCTURE) == 1.0 ){
                    boundary_face_node++;
                    const double contact_angle_i = r_face[i].FastGetSolutionStepValue(CONTACT_ANGLE);
                    if (contact_angle_i > 1.0e-12)
                    {
                        contact_angle += contact_angle_i;
                        contact_angle_weight += 1.0;
                    }
                    solid_normal += r_face[i].FastGetSolutionStepValue(NORMAL);
                }
            }

            if (boundary_face_node == num_face_nodes){
                double minus_cos_contact_angle = 0.0;
                const double norm_solid_normal = Kratos::norm_2(solid_normal);
                solid_normal /= norm_solid_normal;

                const unsigned int num_int_pts = (faces[i_face]).IntegrationPointsNumber(integration_method);

                auto face_gauss_pts = (faces[i_face]).IntegrationPoints(integration_method);

                VectorType face_jacobian;
                (faces[i_face]).DeterminantOfJacobian(face_jacobian, integration_method);

                // Get the original geometry shape function and gradients values over the intersection
                for (unsigned int i_gauss = 0; i_gauss < num_int_pts; ++i_gauss) {

                    // Store the Gauss points weights
                    const double face_weight = face_jacobian(i_gauss) * face_gauss_pts[i_gauss].Weight();

                    // Compute the global coordinates of the face Gauss pt.
                    GeometryType::CoordinatesArrayType global_coords = ZeroVector(num_dim);
                    global_coords = (faces[i_face]).GlobalCoordinates(global_coords, face_gauss_pts[i_gauss].Coordinates());

                    // Compute the parent geometry local coordinates of the Gauss pt.
                    GeometryType::CoordinatesArrayType loc_coords = ZeroVector(num_dim);
                    loc_coords = GetGeometry().PointLocalCoordinates(loc_coords, global_coords);

                    // Compute shape function values
                    // Obtain the parent subgeometry shape function values

                    Kratos::Vector face_shape_func = ZeroVector(num_nodes);
                    face_shape_func = GetGeometry().ShapeFunctionsValues(face_shape_func, loc_coords);

                    for (unsigned int i_node = 0; i_node < num_nodes; i_node++){

                        // Working with averaged unmodified DISTANE_GRADIENT in the 0th-step?
                        //VectorType grad_phi_old_avg_i = GetGeometry()[i_node].FastGetSolutionStepValue(DISTANCE_GRADIENT);
                        //grad_phi_old_avg_i /= norm_2(grad_phi_old_avg_i); //It will help? or disturb? the back of the droplet (upstream where grad < 1) if used without prior parallel redistancing

                        // Working with averaged DISTANE_GRADIENT (accessed in the process) for the 0th-step?
                        VectorType grad_phi_avg_i = GetGeometry()[i_node].FastGetSolutionStepValue(DISTANCE_GRADIENT);//.GetValue(DISTANCE_GRADIENT);
                        const double norm_grad_phi_avg_i = norm_2( grad_phi_avg_i );
                        const double distance_i = GetGeometry()[i_node].FastGetSolutionStepValue(DISTANCE);
                        //grad_phi_avg_i /= norm_2(grad_phi_avg_i); // It is not a good idea!

                        /* const double normalDist = 1.0/(5.0*element_size)* std::abs(GetGeometry()[i_node].GetSolutionStepValue(DISTANCE));
                        double theta = 1.0;
                        if (normalDist < 1.0){
                            theta = 1.0 - std::max(0.0, std::min( normalDist*normalDist*(3.0 - 2.0*normalDist), 1.0 ));//0.0; for small unpinned hydrophobic droplet //1.0; for pinned droplet
                        } */
                        const double theta = 1.0; //std::min( std::exp( -( (std::abs(distance_i) + tolerance)/(20.0*element_size) - 1.0 ) ), 1.0); //1.0;

                        if (contact_angle_weight > 0.0){
                            minus_cos_contact_angle = -theta/* *norm_grad_phi_avg_i */*std::cos(contact_angle/contact_angle_weight) +
                            (1.0-theta)*Kratos::inner_prod(solid_normal,grad_phi_avg_i)/norm_grad_phi_avg_i;
                            /* minus_cos_contact_angle = -std::cos(contact_angle/contact_angle_weight);
                            minus_cos_contact_angle = minus_cos_contact_angle*norm_grad_phi_avg_i; */
                        } else{
                            minus_cos_contact_angle = Kratos::inner_prod(solid_normal,grad_phi_avg_i)/norm_grad_phi_avg_i;
                        }

                        /* if (step == 0){ */
                            rhs(i_node) += scale * minus_cos_contact_angle * face_weight * face_shape_func(i_node);
                        /* } else{ // step == 1
                            rhs(i_node) += minus_cos_contact_angle * face_weight * face_shape_func(i_node);
                        } */
                        //Note: the elliptic redistancing doesn't need BC.
                    }
                }
            }
            i_face++;
        }
    }

    noalias(rLeftHandSideMatrix) = lhs + lhs_Newton_Raphson;
    noalias(rRightHandSideVector) = rhs - prod(lhs,values); // Reducing the contribution of the last known values

    KRATOS_CATCH("");
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix only
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void EDReinitializationElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void EDReinitializationElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the first derivatives contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void EDReinitializationElement::CalculateFirstDerivativesContributions(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix for the first derivatives constributions
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void EDReinitializationElement::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the first derivatives constributions
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void EDReinitializationElement::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * ELEMENTS inherited from this class must implement this methods
 * if they need to add dynamic element contributions
 * note: second derivatives means the accelerations if the displacements are the dof of the analysis
 * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
 * CalculateSecondDerivativesContributions,
 * CalculateSecondDerivativesLHS, CalculateSecondDerivativesRHS methods are : OPTIONAL
 */


/**
 * this is called during the assembling process in order
 * to calculate the second derivative contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void EDReinitializationElement::CalculateSecondDerivativesContributions(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix for the second derivatives constributions
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void EDReinitializationElement::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the second derivatives constributions
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void EDReinitializationElement::CalculateSecondDerivativesRHS(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental mass matrix
 * @param rMassMatrix: the elemental mass matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void EDReinitializationElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (rMassMatrix.size1() != 0)
        rMassMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental damping matrix
 * @param rDampingMatrix: the elemental damping matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void EDReinitializationElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (rDampingMatrix.size1() != 0)
        rDampingMatrix.resize(0, 0, false);
}

/**
 * This method provides the place to perform checks on the completeness of the input
 * and the compatibility with the problem options as well as the contitutive laws selected
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 * this method is: MANDATORY
 */
int EDReinitializationElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->Id() < 1) <<"EllipticDistanceReinitializationElement found with Id 0 or negative" << std::endl;

    KRATOS_ERROR_IF(this->GetGeometry().Area() <= 0) << "On EllipticDistanceReinitializationElement -> "
        << this->Id() <<  "; Area cannot be less than or equal to 0" << std::endl;

      // Base class checks for positive Jacobian and Id > 0
      int ierr = Element::Check(rCurrentProcessInfo);
      if(ierr != 0) return ierr;

      // Check that all required variables have been registered
    //   KRATOS_CHECK_VARIABLE_KEY(DISTANCE);
    //   KRATOS_CHECK_VARIABLE_KEY(DISTANCE_AUX2);

      unsigned const int number_of_points = GetGeometry().size();
      // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
      for ( unsigned int i = 0; i < number_of_points; i++ )
      {
          auto &rnode = this->GetGeometry()[i];//Node<3UL>
          KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISTANCE_AUX2,rnode)
          KRATOS_CHECK_DOF_IN_NODE(DISTANCE_AUX2,rnode)
      }

      return ierr;

    KRATOS_CATCH("");
}

///@}
///@name Access
///@{


///@}
///@name Inquiry
///@{


///@}
///@name Input and output
///@{

/// Turn back information as a string.

std::string EDReinitializationElement::Info() const {
    std::stringstream buffer;
    buffer << "EllipticDistanceReinitializationElement #" << Id();
    return buffer.str();
}

/// Print information about this object.

void EDReinitializationElement::PrintInfo(std::ostream& rOStream) const {
    rOStream << "EllipticDistanceReinitializationElement #" << Id();
}

/// Print object's data.

void EDReinitializationElement::PrintData(std::ostream& rOStream) const {
    pGetGeometry()->PrintData(rOStream);
}

///@}
///@name Friends
///@{

///@}

///@name Protected static Member Variables
///@{

///@}
///@name Protected member Variables
///@{

///@}
///@name Protected Operators
///@{

///@}
///@name Protected Operations
///@{

///@}
///@name Protected  Access
///@{

///@}
///@name Protected Inquiry
///@{

///@}
///@name Protected LifeCycle
///@{

///@}

///@name Static Member Variables
///@{

///@}
///@name Member Variables
///@{

///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{

///@}
///@name Serialization
///@{

void EDReinitializationElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );

    // List
    // To be completed with the class member list
}

void EDReinitializationElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );

    // List
    // To be completed with the class member list
}

///@}
///@name Private  Access
///@{

///@}
///@name Private Inquiry
///@{

///@}
///@name Un accessible methods
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >> (std::istream& rIStream, EDReinitializationElement& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const EDReinitializationElement& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.
