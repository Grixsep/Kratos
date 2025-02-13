// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Joaquin Irazabal Gonzalez
//

// System includes

// External includes


// Project includes
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

// Application includes
#include "free_field_interface_element.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
FreeFieldInterfaceElement::FreeFieldInterfaceElement( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseInterfaceElement( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

FreeFieldInterfaceElement::FreeFieldInterfaceElement( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : BaseInterfaceElement( NewId, pGeometry, pProperties )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer FreeFieldInterfaceElement::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<FreeFieldInterfaceElement>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer FreeFieldInterfaceElement::Create( IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<FreeFieldInterfaceElement>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

FreeFieldInterfaceElement::~FreeFieldInterfaceElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer FreeFieldInterfaceElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    FreeFieldInterfaceElement::Pointer p_new_elem = Kratos::make_intrusive<FreeFieldInterfaceElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // Currently selected integration methods
    p_new_elem->SetIntegrationMethod(BaseType::mThisIntegrationMethod);

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(BaseType::mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

bool FreeFieldInterfaceElement::UseElementProvidedStrain() const
{
    return true;
}

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldInterfaceElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const auto& r_prop = GetProperties();

    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    noalias(rMassMatrix) = ZeroMatrix( mat_size, mat_size );

    // Checking density
    KRATOS_ERROR_IF_NOT(r_prop.Has(DENSITY)) << "DENSITY has to be provided for the calculation of the MassMatrix!" << std::endl;

    const double density = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
    const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;
    const double width = 1.0; // Width of the element, for the moment it is set to 1.0

    // Calculate the length in the y direction between two nodes
    double length = 0.0;
    if (dimension == 2 && number_of_nodes >= 2) {
        const auto& node1 = r_geom[0];
        const auto& node2 = r_geom[1];

        // Node coordinates
        const double y1 = node1.Y0();
        const double y2 = node2.Y0();

        // Calculate the length in the y direction
        length = std::abs(y2 - y1);
    }

    // Total mass of the element
    const double element_mass = density * length * width * thickness;

    // Assemble the mass matrix
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        for (SizeType j = 0; j < dimension; ++j) {
            const SizeType index = i * dimension + j;
            if (index == 4 || index == 5 || index == 6 || index == 7) {
                rMassMatrix(index, index) = 0.5 * element_mass;
            } else {
                rMassMatrix(index, index) = 0.0;
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldInterfaceElement::CalculateDampingMatrix(
    MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const auto& r_prop = GetProperties();

    const auto& r_geom = GetGeometry();
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    unsigned int mat_size = number_of_nodes * dimension;

    if ( rDampingMatrix.size1() != mat_size )
        rDampingMatrix.resize( mat_size, mat_size, false );

    Kratos::noalias( rDampingMatrix ) = ZeroMatrix( mat_size, mat_size );

    // Get material properties
    double wave_velocity_p = this->GetValue(WAVE_VELOCITY_P);
    double wave_velocity_s = this->GetValue(WAVE_VELOCITY_S);

    // Checking density
    KRATOS_ERROR_IF_NOT(r_prop.Has(DENSITY)) << "DENSITY has to be provided for the calculation of the MassMatrix!" << std::endl;

    const double density = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
    const double thickness = (dimension == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;
    const double width = 1.0; // Width of the element, for the moment it is set to 1.0

    // Calculate the length in the y direction between two nodes
    double area = 0.0;
    if (dimension == 2 && number_of_nodes >= 2) {
        const auto& node1 = r_geom[0];
        const auto& node2 = r_geom[1];

        // Node coordinates
        const double y1 = node1.Y0();
        const double y2 = node2.Y0();

        // Calculate the length in the y direction
        area = std::abs(y2 - y1) * width;
    }

    // Calculate damping matrix elements (Asymmetric matrix)
    // Nielsen, A. H. (2006, May). Absorbing boundary conditions for seismic analysis in ABAQUS. In ABAQUS users’ conference (pp. 359-376).
    rDampingMatrix(0, 0) = +0.5 * density * area * thickness * wave_velocity_p; // c33
    rDampingMatrix(2, 2) = +0.5 * density * area * thickness * wave_velocity_p; // c55
    rDampingMatrix(0, 6) = -0.5 * density * area * thickness * wave_velocity_p; // c31
    rDampingMatrix(2, 4) = -0.5 * density * area * thickness * wave_velocity_p; // c57

    rDampingMatrix(1, 1) = +0.5 * density * area * thickness * wave_velocity_s; // c44
    rDampingMatrix(3, 3) = +0.5 * density * area * thickness * wave_velocity_s; // c66
    rDampingMatrix(1, 7) = -0.5 * density * area * thickness * wave_velocity_s; // c42
    rDampingMatrix(3, 5) = -0.5 * density * area * thickness * wave_velocity_s; // c68

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldInterfaceElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType strain_size = GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    const bool is_rotated = IsElementRotated();

    KinematicVariables this_kinematic_variables(strain_size, dimension, number_of_nodes);
    ConstitutiveVariables this_constitutive_variables(strain_size);

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size )
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size() != mat_size )
            rRightHandSideVector.resize( mat_size, false );

        noalias(rRightHandSideVector) = ZeroVector( mat_size ); //resetting RHS
    }

    // Reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = this->IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters Values(r_geometry,GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, UseElementProvidedStrain());
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    if ( CalculateStiffnessMatrixFlag ) {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    } else {
        ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
    }

    // If strain has to be computed inside of the constitutive law with PK2
    Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

    // Some declarations
    array_1d<double, 3> body_force;
    double int_to_reference_weight;

    // Computing in all integrations points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod());

        // Compute material response
        CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, GetStressMeasure(), is_rotated);

        // Calculating weights for integration on the reference configuration
        int_to_reference_weight = GetIntegrationWeight(integration_points, point_number, this_kinematic_variables.detJ0);

        if ( dimension == 2 && GetProperties().Has( THICKNESS ))
            int_to_reference_weight *= GetProperties()[THICKNESS];

        if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
            // Contributions to stiffness matrix calculated on the reference config
            this->CalculateAndAddKm( rLeftHandSideMatrix, this_kinematic_variables.B, this_constitutive_variables.D, int_to_reference_weight );
        }
    }

    if ( CalculateStiffnessMatrixFlag && CalculateResidualVectorFlag ) { // Calculation of the matrix is required
        // Obtain the displacements
        Vector displacements(mat_size);
        GetValuesVector(displacements);

        Vector stiffness_contribution = prod(rLeftHandSideMatrix, displacements);

        for (SizeType i = 0; i < mat_size; ++i) {
            rRightHandSideVector[i] -= stiffness_contribution[i];
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldInterfaceElement::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod
    )
{
    const GeometryType::IntegrationPointsArrayType& r_integration_points = this->IntegrationPoints(rIntegrationMethod);
    // Shape functions
    if (PointNumber == 0 || PointNumber == 1) {
        rThisKinematicVariables.N[0] =  rThisKinematicVariables.N[1] =  rThisKinematicVariables.N[2] =  rThisKinematicVariables.N[3] =  0.0;
    } else {
        rThisKinematicVariables.N[0] =  0.0;
        rThisKinematicVariables.N[1] =  0.0;
        rThisKinematicVariables.N[2] =  0.5 * ( 1.0 - r_integration_points[PointNumber].Coordinates()[1] );
        rThisKinematicVariables.N[3] =  0.5 * ( 1.0 + r_integration_points[PointNumber].Coordinates()[1] );
    }

    rThisKinematicVariables.detJ0 = CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);

    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;

    // Compute B
    CalculateB( rThisKinematicVariables.B, rThisKinematicVariables.DN_DX, r_integration_points, PointNumber );

    // Compute equivalent F
    GetValuesVector(rThisKinematicVariables.Displacements);
    Vector strain_vector(mConstitutiveLawVector[0]->GetStrainSize());
    noalias(strain_vector) = prod(rThisKinematicVariables.B, rThisKinematicVariables.Displacements);
    ComputeEquivalentF(rThisKinematicVariables.F, strain_vector);
    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldInterfaceElement::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;

    // Displacements vector
    Vector displacements(mat_size);
    GetValuesVector(displacements);

    // Compute strain
    noalias(rThisConstitutiveVariables.StrainVector) = prod(rThisKinematicVariables.B, displacements);

    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); //assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); //assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldInterfaceElement::CalculateB(
    Matrix& rB,
    const Matrix& rDN_DX,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const IndexType PointNumber
    ) const
{
    KRATOS_TRY;

    StructuralMechanicsElementUtilities::CalculateB(*this, rDN_DX, rB);

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldInterfaceElement::ComputeEquivalentF(
    Matrix& rF,
    const Vector& rStrainTensor
    ) const
{
    StructuralMechanicsElementUtilities::ComputeEquivalentF(*this, rStrainTensor, rF);
}

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldInterfaceElement::CalculateAndAddKm(
    MatrixType& rLeftHandSideMatrix,
    const Matrix& B,
    const Matrix& D,
    const double IntegrationWeight
    ) const
{
    noalias( rLeftHandSideMatrix ) += IntegrationWeight * prod( trans( B ), Matrix(prod(D, B)));

    // Get material properties
    const double n_dir = this->GetValue(NORMAL_DIRECTION);

    const double young = 34e9;
    const double poisson = 0.29;
    const double lambda = young * poisson / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    const double mu = young / (2.0 * (1.0 + poisson));

    // Asymmetric matrix:
    // Nielsen, A. H. (2006, May). Absorbing boundary conditions for seismic analysis in ABAQUS. In ABAQUS users’ conference (pp. 359-376).
    rLeftHandSideMatrix(0, 7) = +0.5 * n_dir * lambda; // k32
    rLeftHandSideMatrix(2, 7) = +0.5 * n_dir * lambda; // k52
    rLeftHandSideMatrix(0, 5) = -0.5 * n_dir * lambda; // k38
    rLeftHandSideMatrix(2, 5) = -0.5 * n_dir * lambda; // k58

    rLeftHandSideMatrix(1, 6) = +0.5 * n_dir * mu; // k41
    rLeftHandSideMatrix(3, 6) = +0.5 * n_dir * mu; // k61
    rLeftHandSideMatrix(1, 4) = -0.5 * n_dir * mu; // k47
    rLeftHandSideMatrix(3, 4) = -0.5 * n_dir * mu; // k67
}

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldInterfaceElement::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseInterfaceElement );
}

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldInterfaceElement::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseInterfaceElement );
}

} // Namespace Kratos


