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
#include "custom_conditions/point_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

PointLoadCondition::PointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseLoadCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

PointLoadCondition::PointLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : BaseLoadCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer PointLoadCondition::Create(IndexType NewId,GeometryType::Pointer pGeom,PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<PointLoadCondition>(NewId, pGeom, pProperties);
}

//************************************************************************************
//************************************************************************************

Condition::Pointer PointLoadCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_intrusive<PointLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer PointLoadCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<PointLoadCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

PointLoadCondition::~PointLoadCondition()
{
}

//************************************************************************************
//************************************************************************************

void PointLoadCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    const unsigned int NumberOfNodes = GetGeometry().size();
    const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();

    // Resizing as needed the LHS
    const unsigned int MatSize = NumberOfNodes * Dimension;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
        {
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size( ) != MatSize )
        {
            rRightHandSideVector.resize( MatSize, false );
        }

        // noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
    }

     // Reading integration points and local gradients
    const auto integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(GetGeometry());
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(integration_method);
    const Matrix& rNcontainer = GetGeometry().ShapeFunctionsValues(integration_method);

    // Sizing work matrices
    Vector stress_on_nodes_x = ZeroVector( NumberOfNodes );
    Vector stress_on_nodes_y = ZeroVector( NumberOfNodes );
    Vector stress_on_nodes_z = ZeroVector( NumberOfNodes );

    // Pressure applied to the element itself
    if (Dimension == 2) {
        double stress_on_condition_x = 0.0;
        double stress_on_condition_y = 0.0;
        double stress_on_condition_z = 0.0;
        if( this->Has( POINT_STRESS_X ) ) {
            stress_on_condition_x += this->GetValue( POINT_STRESS_X );
        }
        if( this->Has( POINT_STRESS_Y ) ) {
            stress_on_condition_y += this->GetValue( POINT_STRESS_Y );
        }
        if( this->Has( POINT_STRESS_Z ) ) {
            stress_on_condition_z += this->GetValue( POINT_STRESS_Z );
        }

        for ( IndexType i = 0; i < stress_on_nodes_x.size(); i++ ) {
            stress_on_nodes_x[i] = stress_on_condition_x;
            if( GetGeometry()[i].SolutionStepsDataHas( POINT_STRESS_X ) ) {
                stress_on_nodes_x[i] -= GetGeometry()[i].FastGetSolutionStepValue( POINT_STRESS_X );
            }
        }
        for ( IndexType i = 0; i < stress_on_nodes_y.size(); i++ ) {
            stress_on_nodes_y[i] = stress_on_condition_y;
            if( GetGeometry()[i].SolutionStepsDataHas( POINT_STRESS_Y ) ) {
                stress_on_nodes_y[i] -= GetGeometry()[i].FastGetSolutionStepValue( POINT_STRESS_Y );
            }
        }
        for ( IndexType i = 0; i < stress_on_nodes_z.size(); i++ ) {
            stress_on_nodes_z[i] = stress_on_condition_z;
            if( GetGeometry()[i].SolutionStepsDataHas( POINT_STRESS_Z ) ) {
                stress_on_nodes_z[i] -= GetGeometry()[i].FastGetSolutionStepValue( POINT_STRESS_Z );
            }
        }
    }

    // Vector with a loading applied to the condition
    array_1d<double, 3 > PointLoad = ZeroVector(3);
    if( this->Has( POINT_LOAD ) )
    {
        noalias(PointLoad) = this->GetValue( POINT_LOAD );
    }

    for (unsigned int ii = 0; ii < NumberOfNodes; ++ii)
    {
        const unsigned int base = ii*Dimension;

        if( GetGeometry()[ii].SolutionStepsDataHas( POINT_LOAD ) )
        {
            noalias(PointLoad) += GetGeometry()[ii].FastGetSolutionStepValue( POINT_LOAD );
        }

        for(unsigned int k = 0; k < Dimension; ++k)
        {
            rRightHandSideVector[base + k] -= GetPointLoadIntegrationWeight() * PointLoad[k];
        }
    }

    // Declaring tangent and Jacobian
    array_1d<double, 3> tangent_xi, tangent_eta;
    Matrix J0(Dimension, 1);

    // Iterate over the Gauss points
    for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
        GeometryUtils::JacobianOnInitialConfiguration(GetGeometry(), integration_points[point_number], J0);
        const double detJ0 = MathUtils<double>::GeneralizedDet(J0);
        const double integration_weight = this->GetIntegrationWeight(integration_points, point_number, detJ0);

        // Calculating the pressure on the gauss point
        double gauss_stress_x = 0.0;
        double gauss_stress_y = 0.0;
        double gauss_stress_z = 0.0;
        for ( IndexType ii = 0; ii < NumberOfNodes; ii++ ) {
            gauss_stress_x += rNcontainer( point_number, ii ) * stress_on_nodes_x[ii];
            gauss_stress_y += rNcontainer( point_number, ii ) * stress_on_nodes_y[ii];
            gauss_stress_z += rNcontainer( point_number, ii ) * stress_on_nodes_z[ii];
        }

        // Adding contributions to the residual vector
        if ( CalculateResidualVectorFlag ) {
            if ( gauss_stress_x != 0.0 ) {
                array_1d<double, 3> normal = ZeroVector(3);
                normal[0] = 1.0;

                this->CalculateAndAddPressureForce( rRightHandSideVector, row( rNcontainer, point_number ), normal, gauss_stress_x, integration_weight );
            }
            if ( gauss_stress_y != 0.0 ) {
                array_1d<double, 3> normal = ZeroVector(3);
                normal[1] = 1.0;

                this->CalculateAndAddPressureForce( rRightHandSideVector, row( rNcontainer, point_number ), normal, gauss_stress_y, integration_weight );
            }
            if ( gauss_stress_z != 0.0 ) {
                array_1d<double, 3> normal = ZeroVector(3);
                normal[2] = 1.0;

                this->CalculateAndAddPressureForce( rRightHandSideVector, row( rNcontainer, point_number ), normal, gauss_stress_z, integration_weight );
            }
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void PointLoadCondition::CalculateAndAddPressureForce(
    Vector& rRightHandSideVector,
    const Vector& rN,
    const array_1d<double, 3>& rNormal,
    double Pressure,
    double IntegrationWeight
    ) const
{
    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType block_size = this->GetBlockSize();
    const unsigned int Dimension = GetGeometry().WorkingSpaceDimension();

    for ( IndexType i = 0; i < number_of_nodes; ++i ) {
        const IndexType index = block_size * i;

        const double coeff = Pressure * rN[i] * IntegrationWeight;

        for ( IndexType j = 0; j < Dimension; ++j ) {
            rRightHandSideVector[index + j] -= coeff * rNormal[j];
        }
    }
}

//************************************************************************************
//************************************************************************************

double PointLoadCondition::GetPointLoadIntegrationWeight() const
{
    return 1.0;
}

} // Namespace Kratos


