// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Joaquín Irazábal González
//

// System includes

// External includes

// Project includes
#include "custom_conditions/free_field_point_load_condition.h"
#include "custom_utilities/structural_mechanics_math_utilities.hpp"

namespace Kratos
{
//******************************* CONSTRUCTOR ****************************************
//************************************************************************************

FreeFieldPointLoadCondition::FreeFieldPointLoadCondition(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    )
        : PointLoadCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

FreeFieldPointLoadCondition::FreeFieldPointLoadCondition(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    )
        : PointLoadCondition( NewId, pGeometry, pProperties )
{
}

//********************************* CREATE *******************************************
//************************************************************************************

Condition::Pointer FreeFieldPointLoadCondition::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<FreeFieldPointLoadCondition>( NewId, pGeom, pProperties );
}

//************************************************************************************
//************************************************************************************

Condition::Pointer FreeFieldPointLoadCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<FreeFieldPointLoadCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer FreeFieldPointLoadCondition::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<FreeFieldPointLoadCondition>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

//******************************* DESTRUCTOR *****************************************
//************************************************************************************

FreeFieldPointLoadCondition::~FreeFieldPointLoadCondition()
{
}

//************************************************************************************
//********************************* PROTECTED ****************************************
//************************************************************************************

void FreeFieldPointLoadCondition::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType block_size = this->GetBlockSize();

    // Get problem dimensions
    const SizeType n_dim = rCurrentProcessInfo[DOMAIN_SIZE];

    // Resizing as needed the LHS
    const unsigned int mat_size = number_of_nodes * block_size;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != mat_size )
        {
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );
        }

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size( ) != mat_size )
        {
            rRightHandSideVector.resize( mat_size, false );
        }

        noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
    }

    // Vector with a loading applied to the condition
    array_1d<double, 3> normal = ZeroVector(3);
    double density = 0.0;
    double wave_velocity_p = 0.0;
    double wave_velocity_s = 0.0;

    if( this->Has( DENSITY ) )
    {
        density = this->GetValue( DENSITY );
    }

    if( this->Has( WAVE_VELOCITY_P ) )
    {
        wave_velocity_p = this->GetValue( WAVE_VELOCITY_P );
    }

    if( this->Has( WAVE_VELOCITY_S ) )
    {
        wave_velocity_s = this->GetValue( WAVE_VELOCITY_S );
    }

    // Add free-field contributions
    for (IndexType ii = 0; ii < number_of_nodes; ++ii) {
        const auto& velocity = r_geometry[ii].FastGetSolutionStepValue(VELOCITY);

        // Compute damping forces
        array_1d<double, 3> damping_force = ZeroVector(3);

        damping_force[0] = -density * (wave_velocity_s * velocity[0]);
        damping_force[1] = -density * (wave_velocity_p * velocity[1]);
        damping_force[2] = -density * (wave_velocity_s * velocity[2]);

        // Add damping force to RHS
        const IndexType base = ii * block_size;
        for (IndexType k = 0; k < n_dim; ++k) {
            rRightHandSideVector[base + k] += GetPointLoadIntegrationWeight() * damping_force[k];
        }
    }

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void FreeFieldPointLoadCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointLoadCondition );
}

//************************************************************************************
//************************************************************************************

void FreeFieldPointLoadCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointLoadCondition );
}

} // Namespace Kratos


