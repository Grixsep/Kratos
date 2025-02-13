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

        // Declaring normal and tangential vectors
        array_1d<double, 3> normal = ZeroVector(3);
        array_1d<double, 3> tangent_xi = ZeroVector(3);
        array_1d<double, 3>  tangent_eta = ZeroVector(3);

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
            rRightHandSideVector[base + k] += 0.5 * GetPointLoadIntegrationWeight() * damping_force[k];
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


