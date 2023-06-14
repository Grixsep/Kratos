// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    athira vadakkekkara
//  Collaborators:
//

// Project includes
#include "elastic_local_anisotropic_damage_3d.h"
#include "constitutive_laws_application_variables.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "includes/checks.h"
#include "utilities/math_utils.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ElasticAnisotropicDamage::ElasticAnisotropicDamage()
    : ElasticIsotropic3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

ElasticAnisotropicDamage::ElasticAnisotropicDamage(const ElasticAnisotropicDamage &rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer ElasticAnisotropicDamage::Clone() const
{
    return Kratos::make_shared<ElasticAnisotropicDamage>(ElasticAnisotropicDamage(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

ElasticAnisotropicDamage::~ElasticAnisotropicDamage()
{
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

int ElasticAnisotropicDamage::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    KRATOS_CHECK(rMaterialProperties.Has(POISSON_RATIO));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS_TENSION));
    KRATOS_CHECK(rMaterialProperties.Has(YIELD_STRESS_COMPRESSION));
    KRATOS_CHECK(rMaterialProperties.Has(DAMAGE_MODEL_PARAMETER_BETA1_TENSION));
    KRATOS_CHECK(rMaterialProperties.Has(DAMAGE_MODEL_PARAMETER_BETA2_TENSION));
    KRATOS_CHECK(rMaterialProperties.Has(DAMAGE_MODEL_PARAMETER_BETA1_COMPRESSION));
    KRATOS_CHECK(rMaterialProperties.Has(DAMAGE_MODEL_PARAMETER_BETA2_COMPRESSION));
    return 0;
}

//************************************************************************************
//************************************************************************************

bool ElasticAnisotropicDamage::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == STRAIN_ENERGY){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    } else if(rThisVariable == DAMAGE_VARIABLE){
        return false;
    }

    return false;
}

//************************************************************************************
//************************************************************************************

bool ElasticAnisotropicDamage::Has(const Variable<Vector>& rThisVariable)
{
    if(rThisVariable == INTERNAL_VARIABLES){
        // explicitly returning "false", so the element calls CalculateValue(...)s
        return false;
    } else if(rThisVariable == STRAIN){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    } else if(rThisVariable == DAMAGE_VECTOR){
        return true;
    }

    return false;
}

//************************************************************************************
//************************************************************************************

Vector& ElasticAnisotropicDamage::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValues
    )
{
    KRATOS_TRY
    if(rThisVariable == DAMAGE_VECTOR){
        rValues = mDamageVector;
    }
    return rValues;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValues,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY
    if(rThisVariable == DAMAGE_VECTOR){
        mDamageVector = rValues;
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // const double yield_stress = rMaterialProperties[STRESS_LIMITS](0);
    // const double young_modulus = rMaterialProperties[YOUNG_MODULUS];
    // mStrainVariable = yield_stress / std::sqrt(young_modulus);
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    KRATOS_TRY
    Vector damage_vector;
    this->CalculateStressResponse(rParametersValues, damage_vector);
    mDamageVector = damage_vector;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    KRATOS_TRY
    Vector damage_vector;
    CalculateStressResponse(rParametersValues, damage_vector);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::CalculateStressResponse(
    ConstitutiveLaw::Parameters& rParametersValues,
    Vector& rDamageVector)
{
    KRATOS_TRY
    const Properties& r_material_properties = rParametersValues.GetMaterialProperties();
    Flags& r_constitutive_law_options = rParametersValues.GetOptions();
    Vector& r_strain_vector = rParametersValues.GetStrainVector();
    CalculateValue(rParametersValues, STRAIN, r_strain_vector);

    // If we compute the tangent moduli or the stress
    if( r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_STRESS ) ||
        r_constitutive_law_options.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ))
        {
        Vector& r_stress_vector       = rParametersValues.GetStressVector();
        const Vector& r_strain_vector = rParametersValues.GetStrainVector();
        Matrix& r_constitutive_matrix = rParametersValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(r_constitutive_matrix, rParametersValues);
        noalias(r_stress_vector)      = prod(r_constitutive_matrix, r_strain_vector);
        const double beta1t = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA1_TENSION];
        const double beta2t = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA2_TENSION];
        const double beta1c = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA1_COMPRESSION];
        const double beta2c = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA2_COMPRESSION];
        const double E   = r_material_properties[YOUNG_MODULUS];
        const double fck = r_material_properties[YIELD_STRESS_COMPRESSION];
        const double ft  = r_material_properties[YIELD_STRESS_TENSION];
        Vector damage_vector= ZeroVector(3);
        BoundedVectorType Spr = ZeroVector(3);
        BoundedVectorType k0 = ZeroVector(3);
        BoundedVectorType beta1 = ZeroVector(3);
        BoundedVectorType beta2 = ZeroVector(3);
        BoundedVectorType kappa = ZeroVector(3);
        BoundedVectorType F     = ZeroVector(3);
        BoundedMatrixVoigtType EffStiffnessMatrix = ZeroMatrix(6, 6);
        BoundedMatrix3x6Type dEprdE = ZeroMatrix(3,6);
        BoundedMatrixType dkdEpr = ZeroMatrix(3,3);
        BoundedMatrix6x3Type dHdk   = ZeroMatrix(6,3);

        double k0t, k0c;
        if(r_material_properties.Has(DAMAGE_THRESHOLD_TENSION)==true){
            k0t = r_material_properties[DAMAGE_THRESHOLD_TENSION];
        }else{
            k0t = ft/E;
        }
        if(r_material_properties.Has(DAMAGE_THRESHOLD_COMPRESSION)==true){
            k0c = r_material_properties[DAMAGE_THRESHOLD_COMPRESSION];
        }else{
            k0c = (10./3.) * ft/E;
        }
        BoundedVectorType principal_strains = ZeroVector(3);

        GetEigenValues(principal_strains, STRAIN, r_strain_vector);  //calculate prinicpal strains
        for(SizeType i = 0; i < Dimension; ++i) {
            k0[i] = (principal_strains[i] > eps) ? k0t : k0c;
            beta1[i] = (principal_strains[i] > eps) ? beta1t : beta1c;
            beta2[i] = (principal_strains[i] > eps) ? beta2t : beta2c;
            kappa[i]       = std::max(fabs(principal_strains[i]),k0[i]);
            F[i]           = fabs(principal_strains[i])-kappa[i];
        }
        //Compute damage in principal directions
        for (SizeType i = 0; i < Dimension; ++i) {
            if (kappa[i]>= 0 && kappa[i]<=k0[i]) {
                damage_vector[i]=0.0;
            }else if (kappa[i]> k0[i]){
                const double var1      = pow((k0[i]/kappa[i]),beta1[i]);
                const double var2      = exp(-beta2[i]*((kappa[i]-k0[i])/(k0[i])));
                damage_vector[i] = 1.0 - var1 * var2;
            }
            if(damage_vector[i] < 0.0){
                damage_vector[i] = 0.0;
            }
        }
        KRATOS_WATCH(damage_vector)
        CalculateParameters(EffStiffnessMatrix, dEprdE, dkdEpr, rParametersValues, damage_vector);
        noalias(r_stress_vector)       = prod(EffStiffnessMatrix,r_strain_vector);
        CalculatePartialDerivatives(dHdk, r_material_properties, damage_vector, k0, beta1, beta2, kappa);
        const BoundedMatrix3x6Type b = prod(dkdEpr, dEprdE);
        noalias(r_constitutive_matrix) = EffStiffnessMatrix + prod(dHdk, b);
        KRATOS_WATCH(r_constitutive_matrix)
        rDamageVector = damage_vector;
        }
    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::GetEigenValues(
    BoundedVectorType& Pri_Values,
    const Variable<Vector>& rThisVariable,
    const Vector& VectorForm)
{
    KRATOS_TRY
    BoundedMatrixType MatrixForm = ZeroMatrix(3,3);
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    VectorToTensor(MatrixForm, VectorForm, rThisVariable);
    //prinicpal values, max and min
    MathUtils<double>::GaussSeidelEigenSystem(MatrixForm, EigenVectors, EigenValues);
    Pri_Values[0] = EigenValues(0,0);
    Pri_Values[1] = EigenValues(1,1);
    Pri_Values[2] = EigenValues(2,2);
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

Vector& ElasticAnisotropicDamage::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValues
    )
{
    KRATOS_TRY
    if (rThisVariable == DAMAGE_VECTOR) {
        Vector damage_vector;
        this->CalculateStressResponse( rParametersValues, damage_vector);
        rValues = damage_vector;
    }else{
        ElasticIsotropic3D::CalculateValue(rParametersValues, rThisVariable, rValues);
    }
    return( rValues );
    KRATOS_CATCH("")
}

//************************************************************************************
//***********************************************************************************

void ElasticAnisotropicDamage::GetDamageEffectTensor(
    BoundedMatrixVoigtType& DamageEffectTensor,
    const BoundedVectorType& DamageVector

)
{
    KRATOS_TRY
    const double D1 =  DamageVector[0];
    const double D2 =  DamageVector[1];
    const double D3 =  DamageVector[2];
    DamageEffectTensor(0,0) = pow((1-D1),(-1));
    DamageEffectTensor(1,1) = pow((1-D2),(-1));
    DamageEffectTensor(2,2) = pow((1-D3),(-1));
    DamageEffectTensor(3,3) = pow((1-((D2+D3)*0.5)),(-1));
    DamageEffectTensor(4,4) = pow((1-((D3+D1)*0.5)),(-1));
    DamageEffectTensor(5,5) = pow((1-((D1+D2)*0.5)),(-1));

    KRATOS_CATCH("")
}
//************************************************************************************
//***********************************************************************************

void ElasticAnisotropicDamage::CalculateParameters(
    BoundedMatrixVoigtType& EffStiffnessMatrix,
    BoundedMatrix3x6Type& dEprdE,
    BoundedMatrixType& dkdEpr,
    ConstitutiveLaw::Parameters& rParametersValues,
    const Vector& DamageVector
)
{
    KRATOS_TRY
    Matrix& r_constitutive_matrix = rParametersValues.GetConstitutiveMatrix();
    const Vector& r_strain_vector = rParametersValues.GetStrainVector();
    BoundedMatrixVoigtType M     = ZeroMatrix(6,6);
    BoundedMatrixVoigtType Inv_M = ZeroMatrix(6,6);
    BoundedVectorType principal_strains = ZeroVector(3);
    double det_M;
    CalculateElasticMatrix(r_constitutive_matrix, rParametersValues);
    //KRATOS_WATCH(DamageVector)
    GetDamageEffectTensor(M, DamageVector);
    //KRATOS_WATCH(M)
    MathUtils<double>::InvertMatrix(M, Inv_M, det_M);
    //KRATOS_WATCH(Inv_M)
    const BoundedMatrixVoigtType a = prod(r_constitutive_matrix,trans(Inv_M));
    EffStiffnessMatrix = prod(Inv_M,a);
    GetEigenValues(principal_strains, STRAIN, r_strain_vector);
    CalculateDerivativesofEigenvalues(dEprdE, principal_strains, r_strain_vector, STRAIN);
    for(SizeType i = 0; i < Dimension; ++i){
        if(DamageVector[i] > 0){
            dkdEpr(i,i) = 1;
        }
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::CalculatePartialDerivatives(
    BoundedMatrix6x3Type& dHdk,
    const Properties& rMaterialProperties,
    const Vector& DamageVector,
    const BoundedVectorType& Kappa0,
    const BoundedVectorType& Beta1,
    const BoundedVectorType& Beta2,
    const BoundedVectorType& Kappa
)
{
    KRATOS_TRY
    const double E   = rMaterialProperties[YOUNG_MODULUS];
    const double nu  = rMaterialProperties[POISSON_RATIO];
    const double E_factor =  E/((1 + nu ) * (1- 2 * nu));
    BoundedVectorType dDdkappa = ZeroVector(3);
    for(SizeType i =0; i < Dimension; ++i){
        dDdkappa[i]   =  (1- DamageVector[i]) * (Beta1[i]/Kappa[i] + Beta2[i]/Kappa0[i]);
        dHdk(i,i) = E_factor * (-2 * (1-DamageVector[i]) * (1-nu) * dDdkappa[i]);
    }
    dHdk(3,1) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (DamageVector[1]+DamageVector[2]))-1) * dDdkappa[1] ;
    dHdk(3,2) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (DamageVector[1]+DamageVector[2]))-1) * dDdkappa[2];
    dHdk(4,0) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (DamageVector[0]+DamageVector[2]))-1) * dDdkappa[0];
    dHdk(4,2) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (DamageVector[0]+DamageVector[2]))-1) * dDdkappa[2];
    dHdk(5,0) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (DamageVector[0]+DamageVector[1]))-1) * dDdkappa[0];
    dHdk(5,1) = E_factor * 0.5 * (1-2*nu) * ((0.5 * (DamageVector[0]+DamageVector[1]))-1) * dDdkappa[1];
    KRATOS_CATCH("22")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::VectorToTensor(
    BoundedMatrixType& TensorForm,
    const Vector& VectorForm,
    const Variable<Vector>& rThisVariable
)
{
    KRATOS_TRY
    if(rThisVariable == STRESSES){
        TensorForm(0,1)= TensorForm(1,0)= VectorForm[5];
        TensorForm(0,2)= TensorForm(2,0)= VectorForm[4];
        TensorForm(2,1)= TensorForm(1,2)= VectorForm[3];
    }else if(rThisVariable == STRAIN){
        TensorForm(0,1)= TensorForm(1,0)= 0.5 * VectorForm[5];
        TensorForm(0,2)= TensorForm(2,0)= 0.5 * VectorForm[4];
        TensorForm(2,1)= TensorForm(1,2)= 0.5 * VectorForm[3];
    }
    for(SizeType i = 0; i < Dimension; ++i){
        TensorForm(i,i) = VectorForm[i];
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::CalculateDerivativesofEigenvalues(
     BoundedMatrix3x6Type &DerivativesofEigenvalues,
     const BoundedVectorType &EigenvaluesVector,
     const BoundedVectorVoigtType &Voigtform,
     const Variable<Vector>& rThisVariable
    )
{
    KRATOS_TRY
    BoundedMatrixType Matrixform;
    VectorToTensor(Matrixform, Voigtform, rThisVariable);
    BoundedMatrixType DerivtivesMatrix;
    for(SizeType i = 0; i < Dimension; ++i){
        for(SizeType j = 0; j < Dimension; ++j){
            if(i != j && Matrixform(i,j)< eps){
                Matrixform(i,j) = eps;
            }
        }
    }
    for(SizeType i = 0; i < Dimension; ++i){
        BoundedMatrixType AminusLambdaMatrix = Matrixform - EigenvaluesVector[i] * IdentityMatrix(Dimension, Dimension);
        BoundedMatrixType cofactor_matrix = MathUtils<double>::CofactorMatrix(AminusLambdaMatrix);
        const double trace = cofactor_matrix(0,0) + cofactor_matrix(1,1) + cofactor_matrix(2,2);
        DerivtivesMatrix= (1/trace) * cofactor_matrix;
        DerivativesofEigenvalues(i,0) = DerivtivesMatrix(0,0);
        DerivativesofEigenvalues(i,1) = DerivtivesMatrix(1,1);
        DerivativesofEigenvalues(i,2) = DerivtivesMatrix(2,2);
        DerivativesofEigenvalues(i,3) = DerivtivesMatrix(1,2);
        DerivativesofEigenvalues(i,4) = DerivtivesMatrix(0,2);
        DerivativesofEigenvalues(i,5) = DerivtivesMatrix(0,1);
    }
    KRATOS_CATCH("22")
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mOptions.Set(THREE_DIMENSIONAL_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainSize = this->GetStrainSize();
    rFeatures.mSpaceDimension = this->WorkingSpaceDimension();
}

//************************************************************************************
//************************************************************************************



void ElasticAnisotropicDamage::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mDamageVector", mDamageVector);
}

//************************************************************************************
//************************************************************************************

void ElasticAnisotropicDamage::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mDamageVector", mDamageVector);
}



} /* namespace Kratos.*/
