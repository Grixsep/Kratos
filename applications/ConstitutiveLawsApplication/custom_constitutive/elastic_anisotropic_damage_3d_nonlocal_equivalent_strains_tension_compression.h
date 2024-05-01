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
//  Main authors:   athira vadakkekkara
//  Collaborator:

#pragma once

// System includes

// External includes

// Project includes
#include "custom_constitutive/elastic_isotropic_3d.h"

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
 * @class ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC
 * @ingroup StructuralMechanicsApplication
 * @brief
 * @details This material law is defined by the parameters:
 * - YOUNG_MODULUS
 * - POISSON_RATIO
 * - Damage threshold value
 * - Model parameters - Beta1 and Beta2 in Tension and Compression
 * @warning Valid for small strains
 * @note
 * @author Athira
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC
    : public ElasticIsotropic3D
{
public:

    ///@name Type Definitions
    ///@{
    typedef ProcessInfo ProcessInfoType;
    typedef ConstitutiveLaw BaseType;
    typedef std::size_t SizeType;
    /// Static definition of the dimension
    static constexpr SizeType Dimension = 3;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 6;

    // The definition of the bounded vector type
    typedef BoundedVector<double, Dimension > BoundedVectorType;

    // The definition of the bounded vector type
    typedef BoundedVector<double, VoigtSize > BoundedVectorVoigtType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, Dimension, Dimension> BoundedMatrixType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, VoigtSize, VoigtSize> BoundedMatrixVoigtType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, VoigtSize, Dimension> BoundedMatrix6x3Type;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, Dimension, VoigtSize> BoundedMatrix3x6Type;

    struct InternalDamageVariables{
        BoundedMatrixType DamageMatrix = ZeroMatrix(3,3);
        BoundedVectorType StrainHistoryParameter = ZeroVector(3);
        BoundedVectorType kappa = ZeroVector(3);
        Vector EquivalentStrainsTC = ZeroVector(2);
    };

    struct ModelParameters
    {
        Vector Beta1;
        Vector Beta2;
        Vector kappa0;
        Vector kappa;
        Vector strain_history_parameter;
        Vector del_kappa;
        /**
         * The default constructor
         * @param Dimension
         */
        ModelParameters(const SizeType Dimension)
        {
            if (Beta1.size() != Dimension)
                Beta1.resize(Dimension);
            if (Beta2.size() != Dimension)
                Beta2.resize(Dimension);
            if (kappa0.size() != Dimension)
                kappa0.resize(Dimension);
            if (kappa.size() != Dimension)
                kappa.resize(Dimension);
            if (strain_history_parameter.size() != Dimension)
                strain_history_parameter.resize(Dimension);
            if (del_kappa.size() != Dimension)
                del_kappa.resize(Dimension);

            noalias(Beta1)      = ZeroVector(Dimension);
            noalias(Beta2)      = ZeroVector(Dimension);
            noalias(kappa0)     = ZeroVector(Dimension);
            noalias(kappa)      = ZeroVector(Dimension);
            noalias(strain_history_parameter)  = ZeroVector(Dimension);
            noalias(del_kappa)  = ZeroVector(Dimension);

        }
    };
    // Counted pointer
    KRATOS_CLASS_POINTER_DEFINITION(ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC);

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC();

    /**
     * @brief Copy constructor.
     */
    ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC(const ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC& rOther);

    /**
     * @brief Destructor.
     */
    ~ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC() override;

    /**
     * @brief Clone function
     * @return A pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

    ///@}
    ///@name Operators
    ///@{
    ///@}

    ///@name Operations
    ///@{

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    //bool Has(const Variable<double>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    //bool Has(const Variable<Vector>& rThisVariable) override;

    /**
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& GetValue(
        const Variable<Vector>& rThisVariable,
        Vector& rValues
        ) override;


    /**
     * @brief Sets the value of a specified variable (Vector)
     * @param rThisVariable The variable to be returned
     * @param rValue New value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<Vector>& rThisVariable,
        const Vector& rValues,
        const ProcessInfo& rProcessInfo
        ) override;

    /**
     * @brief This is to be called at the very beginning of the calculation
     * @details (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(const Properties& rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector& rShapeFunctionsValues) override;

    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stress
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(Parameters& rValues) override;

    /**
     * @brief Indicates if this CL requires initialization of the material response,
     * called by the element in InitializeMaterialResponse.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return false;
    }

    /**
     * @brief Indicates if this CL requires finalization of the material response,
     * called by the element in FinalizeMaterialResponse.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return true;
    }

    /**
     * @brief Finalize the material response in terms of Cauchy stresses
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;

    // /**
    //  * @brief calculates the value of a specified variable (double)
    //  * @param rValues the needed parameters for the CL calculation
    //  * @param rThisVariable the variable to be returned
    //  * @param rValue a reference to the returned value
    //  * @return rValue output: the value of the specified variable
    //  */
    // double& CalculateValue(Parameters& rValues,
    //                        const Variable<double>& rThisVariable,
    //                        double& rValue) override;

    /**
     * @brief calculates the value of a specified variable (Vector)
     * @param rValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& CalculateValue(Parameters& rValues,
                            const Variable<Vector>& rThisVariable,
                            Vector& rValue) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        rOStream << "Small Strain Isotropic Damage 3D constitutive law\n";
    };

protected:
    ///@name type
    ///@{

    ///@}

    ///@name Protected static Member Variables
    ///@{
    const double eps = 1e-8;
    ///@}

    ///@name Protected member Variables
    ///@{
    Vector mMaterialProperties;
    InternalDamageVariables mInternalDamageVariables;

    ///@}

    ///@name Protected Operators
    ///@{
    ///@}

    ///@name Protected Operations
    ///@{

    /**
     * @brief This method computes principal values of stresses/strains
     * @param VectorForm Stresses/Strains in vector form
     * @param Pri_Values principal values in vector form
     * @param MaxValue maximum of the principal values
     * @param MinValue minimum of the principal values
     */
    void GetEigenValues(BoundedVectorType& Pri_Values,
                        const Variable<Vector>& rThisVariable,
                        const Vector& VectorForm);


    /**
     * @brief This method computes stress weight factor
     */
    void GetStressWeightFactor(double &w,
                               const BoundedVectorType &s_pr
                               ) const ;


    /**
     * @brief This method converts stress or strain vectors to tensors
     */
    void VectorToTensor(BoundedMatrixType& TensorForm,
                        const Vector& VectorForm,
                        const Variable<Vector>& rThisVariable
                        );


    /**
     * @brief calculates the damage effect tensor M
     * @return damage effect tensor M
     */
    void GetDamageEffectTensor(BoundedMatrixVoigtType& DamageEffectTensor,
                               const BoundedVectorType& DamageVector
                               );

    /**
     * @brief Get the Transformed Damageeffect Tensor object
     * @param TransformedDamageEffectTensor
     * @param DamageEffectTensor
     */
    void GetTransformedDamageEffectTensor(BoundedMatrixVoigtType& TransformedDamageEffectTensor,
                                const BoundedMatrixVoigtType& DamageEffectTensor,
                                const Vector& StrainVector
                                );

    /**
     * @brief this method scales the nonlocal equivalent strains to principal direction componenets based on principal strains
     * @param Principal_Nonlocal_Equivalent_Strain
     * @param Principal_Strains
     * @param Nonlocal_Equivalent_Strain
     * @param Local_Equivalent_Strain
     */
    void ScaleNonlocalEquivalentStrain(BoundedVectorType& Principal_Nonlocal_Strains,
                                        ConstitutiveLaw::Parameters& rParametersValues,
                                        const BoundedVectorType& Principal_Strains,
                                        const Vector& Nonlocal_Equivalent_Strains,
                                        const Vector& Local_Equivalent_Strains
                                       );


    /**
     * @brief
     *
     */
    void CalculateLocalEquivalentStrains(Vector& LocalEquivalentStrains,
                                            ConstitutiveLaw::Parameters& rParametersValues,
                                            const BoundedVectorType& PrincipalStrains
                                            );
    /**
     * @brief
     *
     */
    void GetPrincipalDeviatoricStrains(BoundedVectorType& PrincipalDeviatoricStrains,
                                         const Vector& StrainVector,
                                         const BoundedVectorType& PrincipalStrains
                                         );
    /**
     * @brief
     *
     * @param rModelParamaters
     * @param rParametersValues
     * @param PrincipalStrains
     */
    void SetModelParameters(ModelParameters& rModelParamaters,
                            ConstitutiveLaw::Parameters& rParametersValues,
                            const BoundedVectorType& PrincipalStrains,
                            const BoundedVectorType& PrincipalNonlocalStrains
                            );
    void TransformPrincipalDamageIncrementToGlobal(BoundedMatrixType& DamageIncrementTensor,
                                                    const BoundedVectorType& PrinicipalDamageIncrementVector,
                                                    const Vector& StrainVector
                                                    );

    void GetTotalPrincipalDamageVector(BoundedVectorType& PrincipaldamageVector,
                                        const BoundedMatrixType& TotalDamageTensor
                                        );
    /**
     * @brief This method evaluates the Macaulay brackets
     */
    double MacaulayBrackets(const double Number)
    {
        return (Number > 0.0) ? Number : 0.0;
    }

    void ComputePrincipalDamageIncrement(BoundedVectorType& PrincipalDamageIncrement,
                                       ConstitutiveLaw::Parameters& rParametersValues,
                                       ModelParameters rModelParameters,
                                       const BoundedVectorType& PrincipalStrains
                                       );
    void CheckDamageCriteria(
    BoundedVectorType& PrincipalDamageIncrement,
    const ModelParameters& rModelParamaters,
    const BoundedVectorType& PrincipalDamageVector_n
    );
    /**
     * @brief
     * @param H_uNL
     * @param rParametersValues
     * @param Damage_Vector
     * @param Kappa
     * @param local_equivalent_strain
     */
    void Calculate_tangent_HuNL(BoundedVectorVoigtType& H_uNL1,
                                BoundedVectorVoigtType& H_uNL2,
                                ConstitutiveLaw::Parameters& rParametersValues,
                                const ModelParameters& rModelParameters,
                                const Vector& Damage_Vector,
                                const BoundedVectorType& Principal_Strains
                                );


    /**
     * @brief
     * @param H_NLu
     * @param rParametersValues
     * @param Principal_Strains
     */
    void Calculate_tangent_HNLu(BoundedVectorVoigtType& H_NL1u,
                                BoundedVectorVoigtType& H_NL2u,
                                ConstitutiveLaw::Parameters& rParametersValues,
                                const BoundedVectorType& Principal_Strains
                                );

    /**
     * @brief the derivatives of eigen values with respect to the matrix elements
     * @param DerivativesofEigenvalues
     * @param EigenvaluesVector
     * @param Voigtform
     */
    void CalculateDerivativesofEigenvalues(BoundedMatrix3x6Type &DerivativesofEigenvalues,
                                           const BoundedVectorType &EigenvaluesVector,
                                           const BoundedVectorVoigtType &Voigtform,
                                           const Variable<Vector>& rThisVariable);



    /**
     * @brief This method assembles constitutive matrix
     * @param ConstitutiveMatrix
     * @param H_uu
     * @param H_NLu
     * @param H_uNL
     * @param H_NLNL
     */
    void AssembleConstitutiveMatrix(Matrix& ConstitutiveMatrix,
                                    const Matrix& H_uu,
                                    const Vector& H_NL1u,
                                    const Vector& H_NL2u,
                                    const Vector& H_uNL1,
                                    const Vector& H_uNL2,
                                    const double& H_NLNL
                                    );


    

    /**
     * @brief This method computes the stress and constitutive tensor
     * @param rValues The norm of the deviation stress
     * @param rStrainVariable
     */
    void CalculateStressResponse(ConstitutiveLaw::Parameters& rParametersValues,
                                 InternalDamageVariables &rInternalDamageVariables) ;

    ///@}
private:

    ///@name Static Member Variables
    ///@{

    ///@}

    ///@name Member Variables
    ///@{
    
    ///@}

    ///@name Private Operators
    ///@}

    ///@name Private  Access
    ///@{
    ///@}

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // class ElasticAnisotropicDamage3DNonLocalEquivalentStrainsTC
} // namespace Kratos
