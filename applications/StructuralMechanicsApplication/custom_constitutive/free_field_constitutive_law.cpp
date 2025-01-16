// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Joaquin Irazabal González
//
// System includes

// External includes

// Project includes
#include "custom_constitutive/free_field_constitutive_law.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

FreeFieldConstitutiveLaw::FreeFieldConstitutiveLaw()
    : ElasticIsotropic3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

FreeFieldConstitutiveLaw::FreeFieldConstitutiveLaw(const FreeFieldConstitutiveLaw& rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer FreeFieldConstitutiveLaw::Clone() const
{
    return Kratos::make_shared<FreeFieldConstitutiveLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

FreeFieldConstitutiveLaw::~FreeFieldConstitutiveLaw()
{
}

/***********************************************************************************/
/***********************************************************************************/

bool& FreeFieldConstitutiveLaw::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& FreeFieldConstitutiveLaw::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& FreeFieldConstitutiveLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

double& FreeFieldConstitutiveLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
    {
        return rValue;
    }

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void FreeFieldConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( PLANE_STRAIN_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = 3;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 2;
    }

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldConstitutiveLaw::CalculateElasticMatrix(VoigtSizeMatrixType& rC, ConstitutiveLaw::Parameters& rValues)
{
    const auto &r_props = rValues.GetMaterialProperties();
    const double E = r_props[YOUNG_MODULUS];
    const double NU = r_props[POISSON_RATIO];
    ConstitutiveLawUtilities<3>::CalculateElasticMatrixPlaneStrain(rC, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldConstitutiveLaw::CalculatePK2Stress(
    const ConstitutiveLaw::StrainVectorType& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    const double c0 = E / ((1.0 + NU)*(1.0 - 2.0 * NU));
    const double c1 = (1.0 - NU)*c0;
    const double c2 = c0 * NU;
    const double c3 = (0.5 - NU)*c0;

    rStressVector[0] = c1 * rStrainVector[0] + c2 * rStrainVector[1];
    rStressVector[1] = c2 * rStrainVector[0] + c1 * rStrainVector[1];
    rStressVector[2] = c3 * rStrainVector[2];
}

/***********************************************************************************/
/***********************************************************************************/

void FreeFieldConstitutiveLaw::CalculateCauchyGreenStrain(Parameters& rValues, ConstitutiveLaw::StrainVectorType& rStrainVector)
{
    ConstitutiveLawUtilities<3>::CalculateCauchyGreenStrain(rValues, rStrainVector);
}

} // Namespace Kratos
