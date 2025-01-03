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

    const auto& r_geometry = GetGeometry();
    SizeType dimension = r_geometry.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geometry.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    noalias(rMassMatrix) = ZeroMatrix( mat_size, mat_size );

    // Mass Matrix is Zero Matrix

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

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;

    // Asegúrate de que la matriz de amortiguamiento tenga el tamaño correcto
    if (rDampingMatrix.size1() != mat_size || rDampingMatrix.size2() != mat_size)
        rDampingMatrix.resize(mat_size, mat_size, false);

    noalias(rDampingMatrix) = ZeroMatrix(mat_size, mat_size);

    // Obtener las propiedades del elemento
    const double density = GetProperties()[DENSITY];
    const double cp = GetProperties()[WAVE_VELOCITY_P];
    const double cs = GetProperties()[WAVE_VELOCITY_S];

    // Calculate element length
    const double length = r_geometry.Length();

    // Calcular los términos específicos de la matriz de amortiguamiento
    rDampingMatrix(2, 2) = 0.5 * length * density * cp; // c22
    rDampingMatrix(4, 4) = 0.5 * length * density * cp; // c44
    rDampingMatrix(2, 0) = -0.5 * length * density * cp; // c20
    rDampingMatrix(4, 6) = -0.5 * length * density * cp; // c46
    rDampingMatrix(3, 3) = 0.5 * length * density * cs; // c33
    rDampingMatrix(5, 5) = 0.5 * length * density * cs; // c55
    rDampingMatrix(3, 1) = -0.5 * length * density * cs; // c31
    rDampingMatrix(4, 4) = -0.5 * length * density * cs; // c57

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

        if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
            this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, int_to_reference_weight);
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
    const auto& r_geometry = GetGeometry();

    const GeometryType::IntegrationPointsArrayType& r_integration_points = this->IntegrationPoints(rIntegrationMethod);
    // Shape functions
    rThisKinematicVariables.N = r_geometry.ShapeFunctionsValues(rThisKinematicVariables.N, r_integration_points[PointNumber].Coordinates());

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
    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;

    // Asegúrate de que la matriz de rigidez tenga el tamaño correcto
    if (rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size)
        rLeftHandSideMatrix.resize(mat_size, mat_size, false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);

    // Asymmetric matrix:
    // Nielsen, A. H. (2006, May). Absorbing boundary conditions for seismic analysis in ABAQUS. In ABAQUS users’ conference (pp. 359-376).
    rLeftHandSideMatrix(0, 0) = IntegrationWeight * B(0, 0) * D(0, 0); // k00
    rLeftHandSideMatrix(0, 6) = IntegrationWeight * B(0, 6) * D(0, 6); // k06
    rLeftHandSideMatrix(1, 1) = IntegrationWeight * B(1, 1) * D(1, 1); // k11
    rLeftHandSideMatrix(1, 7) = IntegrationWeight * B(1, 7) * D(1, 7); // k17
    rLeftHandSideMatrix(2, 1) = IntegrationWeight * B(2, 1) * D(2, 1); // k21
    rLeftHandSideMatrix(2, 7) = IntegrationWeight * B(2, 7) * D(2, 7); // k27
    rLeftHandSideMatrix(3, 0) = IntegrationWeight * B(3, 0) * D(3, 0); // k30
    rLeftHandSideMatrix(3, 6) = IntegrationWeight * B(3, 6) * D(3, 6); // k36
    rLeftHandSideMatrix(4, 1) = IntegrationWeight * B(4, 1) * D(4, 1); // k41
    rLeftHandSideMatrix(4, 7) = IntegrationWeight * B(4, 7) * D(4, 7); // k47
    rLeftHandSideMatrix(5, 0) = IntegrationWeight * B(5, 0) * D(5, 0); // k50
    rLeftHandSideMatrix(5, 6) = IntegrationWeight * B(5, 6) * D(5, 6); // k56
    rLeftHandSideMatrix(6, 0) = IntegrationWeight * B(6, 0) * D(6, 0); // k60
    rLeftHandSideMatrix(6, 6) = IntegrationWeight * B(6, 6) * D(6, 6); // k66
    rLeftHandSideMatrix(7, 1) = IntegrationWeight * B(7, 1) * D(7, 1); // k71
    rLeftHandSideMatrix(7, 7) = IntegrationWeight * B(7, 7) * D(7, 7); // k77
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


