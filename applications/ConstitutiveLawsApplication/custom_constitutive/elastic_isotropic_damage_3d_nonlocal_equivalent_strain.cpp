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
//  Main authors:   Athira Vadakkekkara
//  Collaborators:
//

// Project includes
#include "elastic_isotropic_damage_3d_nonlocal_equivalent_strain.h"
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

ElasticIsotropicDamage3DNonLocalEquivalentStrain::ElasticIsotropicDamage3DNonLocalEquivalentStrain()
    : ElasticIsotropic3D()
{
}

//********************************COPY CONSTRUCTOR************************************
//************************************************************************************

ElasticIsotropicDamage3DNonLocalEquivalentStrain::ElasticIsotropicDamage3DNonLocalEquivalentStrain(const ElasticIsotropicDamage3DNonLocalEquivalentStrain &rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer ElasticIsotropicDamage3DNonLocalEquivalentStrain::Clone() const
{
    return Kratos::make_shared<ElasticIsotropicDamage3DNonLocalEquivalentStrain>(ElasticIsotropicDamage3DNonLocalEquivalentStrain(*this));
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

ElasticIsotropicDamage3DNonLocalEquivalentStrain::~ElasticIsotropicDamage3DNonLocalEquivalentStrain()
{
}

//********************************DESTRUCTOR******************************************
//************************************************************************************

int ElasticIsotropicDamage3DNonLocalEquivalentStrain::Check(
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

bool ElasticIsotropicDamage3DNonLocalEquivalentStrain::Has(const Variable<double>& rThisVariable)
{
    if(rThisVariable == STRAIN_ENERGY){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    } else if(rThisVariable == EQUIVALENT_STRAIN){
        return true;
    }

    return false;
}

//************************************************************************************
//************************************************************************************

bool ElasticIsotropicDamage3DNonLocalEquivalentStrain::Has(const Variable<Vector>& rThisVariable)
{
    if(rThisVariable == INTERNAL_VARIABLES){
        // explicitly returning "false", so the element calls CalculateValue(...)s
        return false;
    } else if(rThisVariable == STRAIN){
        // explicitly returning "false", so the element calls CalculateValue(...)
        return false;
    }

    return false;
}

//************************************************************************************
//************************************************************************************

double& ElasticIsotropicDamage3DNonLocalEquivalentStrain::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    KRATOS_TRY
    if(rThisVariable == EQUIVALENT_STRAIN){
        rValue = mEquivalentStrain;
    }
    return rValue;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropicDamage3DNonLocalEquivalentStrain::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rProcessInfo
    )
{
    KRATOS_TRY
    if(rThisVariable == EQUIVALENT_STRAIN){
        mEquivalentStrain = rValue;
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropicDamage3DNonLocalEquivalentStrain::InitializeMaterial(
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

void ElasticIsotropicDamage3DNonLocalEquivalentStrain::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    KRATOS_TRY
    double equivalent_strain;
    this->CalculateStressResponse(rParametersValues, equivalent_strain);
    mEquivalentStrain = equivalent_strain;

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropicDamage3DNonLocalEquivalentStrain::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rParametersValues)
{
    KRATOS_TRY
    double equivalent_strain;
    CalculateStressResponse(rParametersValues, equivalent_strain);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropicDamage3DNonLocalEquivalentStrain::CalculateStressResponse(
    ConstitutiveLaw::Parameters& rParametersValues,
    double& rEquivalentStrain)
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
        Matrix r_elastic_tensor;
        r_elastic_tensor.resize(6, 6, false);
        CalculateElasticMatrix(r_elastic_tensor, rParametersValues);
        noalias(r_stress_vector)      = prod(r_elastic_tensor, r_strain_vector);
        BoundedMatrix3x6Type derivatives_of_eigen_values;
        BoundedVectorVoigtType H_uNL        = ZeroVector(6);
        BoundedVectorVoigtType H_NLu        = ZeroVector(6);
        BoundedVectorType Spr = ZeroVector(3);
        BoundedVectorType principal_strains = ZeroVector(3);
        double r, D, local_equivalent_strain, k0t, k0c, SprMax, max_principal_strain;
        const auto& N    = rParametersValues.GetShapeFunctionsValues();
        const double H_NLNL   = 0.0;
        const double eps = 1e-8;
        const double E   = r_material_properties[YOUNG_MODULUS];
        const double fck = r_material_properties[YIELD_STRESS_COMPRESSION];
        const double ft  = r_material_properties[YIELD_STRESS_TENSION];
        GetEigenValues(Spr, SprMax, STRESSES, r_stress_vector);
        const double H = (SprMax < eps) ? 0.0 : 1.0;
        GetEigenValues(principal_strains, max_principal_strain, STRAIN, r_strain_vector);

        //nonlocal equivalent strain
        double nonlocal_equivalent_strain = 0.0;
        for (SizeType i = 0; i < N.size(); ++i) {
            nonlocal_equivalent_strain += N[i] * rParametersValues.GetElementGeometry()[i].FastGetSolutionStepValue(NON_LOCAL_VARIABLE, 0);
        }
        //local equivalent strain
        local_equivalent_strain = sqrt(pow(MacaulayBrackets(principal_strains[0]),2) + pow(MacaulayBrackets(principal_strains[1]),2) + pow(MacaulayBrackets(principal_strains[2]),2));
        GetStressWeightFactor(r,Spr);
        const double del_r = (r > 0 && r < 0.1) ? 1.0 : 0.0;
        const double beta1t = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA1_TENSION];
        const double beta2t = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA2_TENSION];
        const double beta1c = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA1_COMPRESSION];
        const double beta2c = r_material_properties[DAMAGE_MODEL_PARAMETER_BETA2_COMPRESSION];
        if(r_material_properties.Has(DAMAGE_THRESHOLD_TENSION)==true){
            k0t = r_material_properties[DAMAGE_THRESHOLD_TENSION];
        }else{
            k0t = ft/E;
        }
        if(r_material_properties.Has(DAMAGE_THRESHOLD_COMPRESSION)==true){
            k0c = r_material_properties[DAMAGE_THRESHOLD_COMPRESSION];
        }else{
            k0c = ft/E;
        }
        const double k0 = k0t  * H * (1-del_r) + (1.0- H + del_r) * k0c;
        const double beta1 = beta1t  * H * (1-del_r) + (1.-H + del_r) * beta1c;
        const double beta2 = beta2t  * H * (1-del_r) + (1.-H + del_r) * beta2c;
        const double kappa = std::max(nonlocal_equivalent_strain,k0);
        //damage loading condition
        const double f_d   = nonlocal_equivalent_strain - kappa;
        if (f_d < 0.0){
            D = 0.0;
            AssembleConstitutiveMatrix(r_constitutive_matrix, r_elastic_tensor, H_NLu, H_uNL, H_NLNL);
        }else if (f_d == 0.0) {
            const double var1     = pow((k0/kappa),beta1);
            const double var2     = exp(-beta2*((kappa-k0)/(k0)));
            D                     = 1.0 - var1 * var2;
            const double DN       = pow((1.0 - D),2);
            r_stress_vector      *= DN;
            r_elastic_tensor     *= DN ;
            //H_uu = dSigmadEps; H_uNL = dSigmadNL; H_NLu = dlocaldEps; H_NLNL = dlocal/dNL;
            const double dDdKappa = (1.- D) * (beta1/kappa + beta2/k0);
            H_uNL                 = 2.0*(D-1.0) * prod(r_elastic_tensor, r_strain_vector) * dDdKappa;
            CalculateDerivativesofEigenvalues(derivatives_of_eigen_values, principal_strains, r_strain_vector, STRAIN);
            for(SizeType i = 0; i < Dimension; ++i){
                if (MacaulayBrackets(principal_strains[i])!=0.0){
                    for (SizeType j = 0; j < VoigtSize; ++j){
                        H_NLu[j]   +=  principal_strains[i] * derivatives_of_eigen_values(i,j);
                    }
                }
            }
            H_NLu                 /= local_equivalent_strain ;
            AssembleConstitutiveMatrix(r_constitutive_matrix, r_elastic_tensor, H_NLu, H_uNL, H_NLNL);
        }else{
            KRATOS_ERROR << "check the damage loading function" << std::endl;
        }
        if(D < 0.0) D = 0.0;
        rEquivalentStrain = local_equivalent_strain;
    }
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropicDamage3DNonLocalEquivalentStrain::GetEigenValues(
    BoundedVectorType& Principal_Values_vector,
    double& MaxValue,
    const Variable<Vector>& rThisVariable,
    const Vector& VectorForm)
{
    KRATOS_TRY
    BoundedMatrixType MatrixForm = ZeroMatrix(3,3);
    BoundedMatrixType EigenVectors;
    BoundedMatrixType EigenValues = ZeroMatrix(3,3);
    VectorToTensor(MatrixForm, VectorForm, rThisVariable);
    MathUtils<double>::GaussSeidelEigenSystem(MatrixForm, EigenVectors, EigenValues);
    Principal_Values_vector[0] = EigenValues(0,0);
    Principal_Values_vector[1] = EigenValues(1,1);
    Principal_Values_vector[2] = EigenValues(2,2);
    MaxValue    = std::max({Principal_Values_vector[0],Principal_Values_vector[1],Principal_Values_vector[2]}, [](const double a, const double b){return std::abs(b) > std::abs(a);});
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropicDamage3DNonLocalEquivalentStrain::AssembleConstitutiveMatrix(
    Matrix& ConstitutiveMatrix,
    const Matrix& H_uu,
    const Vector& H_NLu,
    const Vector& H_uNL,
    const double& H_NLNL
    )const
{
    KRATOS_TRY
    const SizeType num_rows_C = H_uu.size1();
    const SizeType num_cols_C = H_uu.size2();
    const SizeType size_CuNL  = H_uNL.size();
    const SizeType size_CNLu  = H_NLu.size();

    for (SizeType i = 0; i < num_rows_C ; ++i) {
        for (SizeType j = 0; j < num_cols_C ; ++j) {
            ConstitutiveMatrix(i, j) = H_uu(i,j);
        }
    }
    for (SizeType i = 0; i < size_CuNL; ++i) {
        ConstitutiveMatrix(i, num_cols_C) = H_uNL[i] ;
    }
    for (SizeType i = 0; i < size_CNLu; ++i) {
        ConstitutiveMatrix(num_rows_C, i) = H_NLu[i] ;
    }
    ConstitutiveMatrix(num_rows_C, num_cols_C) = H_NLNL;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropicDamage3DNonLocalEquivalentStrain::VectorToTensor(
    BoundedMatrixType& TensorForm,
    const Vector& VectorForm,
    const Variable<Vector>& rThisVariable
)
{
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
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropicDamage3DNonLocalEquivalentStrain::GetStressWeightFactor(
    double &w,
    const BoundedVectorType &s_pr) const
{
    Vector N1(3);
    const double eps = 1.e-8;
    SizeType kk;
    for(kk=0; kk < Dimension; ++kk ){
        N1(kk) = 0.5 * ( abs(s_pr(kk)) + s_pr(kk) );
    }
    double N11 = N1(0) + N1(1) + N1(2);
    double D11 = eps + abs(s_pr(0)) + abs(s_pr(1)) + abs(s_pr(2));
    w = N11 / D11 ;
}
//************************************************************************************
//************************************************************************************

void ElasticIsotropicDamage3DNonLocalEquivalentStrain::CalculateDerivativesofEigenvalues(
     BoundedMatrix3x6Type &DerivativesofEigenvalues,
     const BoundedVectorType &EigenvaluesVector,
     const BoundedVectorVoigtType &Voigtform,
     const Variable<Vector>& rThisVariable
    )
{
    BoundedMatrixType Matrixform;
    const double eps = 1e-8;
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
}


//************************************************************************************
//************************************************************************************

double& ElasticIsotropicDamage3DNonLocalEquivalentStrain::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    KRATOS_TRY
    if (rThisVariable == EQUIVALENT_STRAIN) {
        double equivalent_strain;
        this->CalculateStressResponse( rParametersValues, equivalent_strain);
        rValue = equivalent_strain;
    }
    return( rValue );
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

Vector& ElasticIsotropicDamage3DNonLocalEquivalentStrain::CalculateValue(
    ConstitutiveLaw::Parameters& rParametersValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    KRATOS_TRY
    //Explicitly having STRAIN and INITAL_STRAIN_VECTOR calculated in base class
    ElasticIsotropic3D::CalculateValue(rParametersValues, rThisVariable, rValue);
    return(rValue);
    KRATOS_CATCH("")
}

//************************************************************************************
//***********************************************************************************

void ElasticIsotropicDamage3DNonLocalEquivalentStrain::GetLawFeatures(Features& rFeatures)
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



void ElasticIsotropicDamage3DNonLocalEquivalentStrain::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.save("mEquivalentStrain", mEquivalentStrain);
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropicDamage3DNonLocalEquivalentStrain::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw);
    rSerializer.load("mEquivalentStrain", mEquivalentStrain);
}



} /* namespace Kratos.*/
