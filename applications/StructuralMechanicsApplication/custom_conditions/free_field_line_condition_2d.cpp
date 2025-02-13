// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_conditions/free_field_line_condition_2d.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "utilities/beam_math_utilities.hpp"
#include "utilities/integration_utilities.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"

namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

FreeFieldLineCondition2D::FreeFieldLineCondition2D(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    ) : LineLoadCondition<2>( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

FreeFieldLineCondition2D::FreeFieldLineCondition2D(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    ) : LineLoadCondition<2>( NewId, pGeometry, pProperties )
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

Condition::Pointer FreeFieldLineCondition2D::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<FreeFieldLineCondition2D>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer FreeFieldLineCondition2D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<FreeFieldLineCondition2D>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer FreeFieldLineCondition2D::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<FreeFieldLineCondition2D>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

/******************************* DESTRUCTOR ****************************************/
/***********************************************************************************/

FreeFieldLineCondition2D::~FreeFieldLineCondition2D()
{
}

/***********************************************************************************/
/********************************* PROTECTED ***************************************/
/***********************************************************************************/

void FreeFieldLineCondition2D::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType block_size = this->GetBlockSize();

    // Get problem dimensions
    const SizeType n_dim = rCurrentProcessInfo[DOMAIN_SIZE];

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * block_size;

    if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size ) {
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size( ) != mat_size ) {
            rRightHandSideVector.resize( mat_size, false );
        }
        noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
    }


    // Reading integration points and local gradients
    const IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geometry);
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(integration_method);

    // Declaring tangent and Jacobian
    array_1d<double, 3> tangent_xi, tangent_eta;
    Matrix J(n_dim, 1);

    // Iterate over the Gauss points
    for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
        const double det_j = r_geometry.DeterminantOfJacobian(integration_points[point_number]);
        const double integration_weight = GetIntegrationWeight(integration_points, point_number, det_j);

        // Compute tangents and normal
        r_geometry.Jacobian(J, point_number, integration_method);
        GetLocalAxis1(tangent_xi, J);
        GetLocalAxis2(tangent_eta);

        tangent_xi /= MathUtils<double>::Norm(tangent_xi);
        tangent_eta /= MathUtils<double>::Norm(tangent_eta);

        array_1d<double, 3> normal;
        MathUtils<double>::UnitCrossProduct(normal, tangent_xi, tangent_eta);

        double density = this->GetValue(DENSITY);
        double wave_velocity_p = this->GetValue(WAVE_VELOCITY_P);
        double wave_velocity_s = this->GetValue(WAVE_VELOCITY_S);

        // Add free-field contributions
        for (IndexType ii = 0; ii < number_of_nodes; ++ii) {
            const auto& velocity = r_geometry[ii].FastGetSolutionStepValue(VELOCITY);

            if (MathUtils<double>::Norm(velocity) < 1e-12) {
                continue; // Skip nodes with zero velocity
            }

            normal[0] = 0.0; normal[1] = 1.0; normal[2] = 0.0;
            tangent_xi[0] = 1.0; tangent_xi[1] = 0.0; tangent_xi[2] = 0.0;
            tangent_eta[0] = 0.0; tangent_eta[1] = 0.0; tangent_eta[2] = 1.0;

            // Compute normal and tangential components
            const double vn = MathUtils<double>::Dot(velocity, normal);
            const double vt1 = MathUtils<double>::Dot(velocity, tangent_xi);
            const double vt2 = MathUtils<double>::Dot(velocity, tangent_eta);

            // Compute damping forces
            array_1d<double, 3> damping_force = -density * (wave_velocity_p * vn * normal + wave_velocity_s * vt1 * tangent_xi + wave_velocity_s * vt2 * tangent_eta);

            // Add damping force to RHS
            const IndexType base = ii * block_size;
            for (IndexType k = 0; k < n_dim; ++k) {
                rRightHandSideVector[base + k] += integration_weight * damping_force[k];
            }
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldLineCondition2D::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LineLoadCondition<2> );
}

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldLineCondition2D::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LineLoadCondition<2> );
}

} // Namespace Kratos




