// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

#pragma once

#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

#include "generalized_newmark_scheme.hpp"

namespace Kratos
{

template <class TSparseSpace, class TDenseSpace>
class NewmarkQuasistaticUPwScheme : public GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(NewmarkQuasistaticUPwScheme);
    using BaseType              = GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>;
    using DofsArrayType         = typename BaseType::DofsArrayType;
    using TSystemMatrixType     = typename BaseType::TSystemMatrixType;
    using TSystemVectorType     = typename BaseType::TSystemVectorType;
    using LocalSystemVectorType = typename BaseType::LocalSystemVectorType;
    using LocalSystemMatrixType = typename BaseType::LocalSystemMatrixType;

    NewmarkQuasistaticUPwScheme(double beta, double gamma, double theta)
        : GeneralizedNewmarkScheme<TSparseSpace, TDenseSpace>(
              {FirstOrderScalarVariable(WATER_PRESSURE, DT_WATER_PRESSURE, DT_PRESSURE_COEFFICIENT)},
              {SecondOrderVectorVariable(DISPLACEMENT), SecondOrderVectorVariable(ROTATION)},
              beta,
              gamma,
              theta)
    {
    }

    void Predict(ModelPart& rModelPart, DofsArrayType& rDofSet, TSystemMatrixType& A, TSystemVectorType& Dx, TSystemVectorType& b) override
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            for (const auto& r_first_order_scalar_variable : this->GetFirstOrderScalarVariables()) {
                if (!rNode.SolutionStepsDataHas(r_first_order_scalar_variable.instance)) continue;
                const auto& instance = r_first_order_scalar_variable.instance;

                if (!rNode.HasDofFor(instance)) continue;

                const auto& first_time_derivative = r_first_order_scalar_variable.first_time_derivative;

                const double previous_variable = rNode.FastGetSolutionStepValue(instance, 1);
                const double current_first_time_derivative =
                    rNode.FastGetSolutionStepValue(first_time_derivative, 0);
                const double previous_first_time_derivative =
                    rNode.FastGetSolutionStepValue(first_time_derivative, 1);
                if (rNode.IsFixed(first_time_derivative)) {
                    rNode.FastGetSolutionStepValue(instance) =
                        previous_variable +
                        this->GetDeltaTime() * this->GetTheta() * previous_first_time_derivative;
                } else if (!rNode.IsFixed(instance)) {
                    rNode.FastGetSolutionStepValue(instance) =
                        previous_variable + this->GetDeltaTime() * previous_first_time_derivative;
                }
            }
        });

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            const std::vector<std::string> components = {"X", "Y", "Z"};
            for (const auto& r_second_order_vector_variable : this->GetSecondOrderVectorVariables()) {
                if (!rNode.SolutionStepsDataHas(r_second_order_vector_variable.instance)) continue;
                for (const auto& component : components) {
                    const auto& instance_component = VariablesUtilities::GetComponentFromVectorVariable(
                        r_second_order_vector_variable.instance.Name(), component);

                    if (!rNode.HasDofFor(instance_component)) continue;

                    const auto& first_time_derivative_component = VariablesUtilities::GetComponentFromVectorVariable(
                        r_second_order_vector_variable.first_time_derivative.Name(), component);
                    const auto& second_time_derivative_component = VariablesUtilities::GetComponentFromVectorVariable(
                        r_second_order_vector_variable.second_time_derivative.Name(), component);

                    const double previous_variable = rNode.FastGetSolutionStepValue(instance_component, 1);
                    const double current_first_time_derivative =
                        rNode.FastGetSolutionStepValue(first_time_derivative_component, 0);
                    const double previous_first_time_derivative =
                        rNode.FastGetSolutionStepValue(first_time_derivative_component, 1);
                    const double current_second_time_derivative =
                        rNode.FastGetSolutionStepValue(second_time_derivative_component, 0);
                    const double previous_second_time_derivative =
                        rNode.FastGetSolutionStepValue(second_time_derivative_component, 1);
                    if (rNode.IsFixed(second_time_derivative_component)) {
                        rNode.FastGetSolutionStepValue(instance_component) =
                            previous_variable + this->GetDeltaTime() * previous_first_time_derivative +
                            this->GetDeltaTime() * this->GetDeltaTime() *
                                ((0.5 - this->GetBeta()) * previous_second_time_derivative +
                                 this->GetBeta() * current_second_time_derivative);
                    } else if (rNode.IsFixed(first_time_derivative_component)) {
                        rNode.FastGetSolutionStepValue(instance_component) =
                            previous_variable +
                            this->GetDeltaTime() *
                                ((this->GetBeta() / this->GetGamma()) *
                                     (current_first_time_derivative - previous_first_time_derivative) +
                                 previous_first_time_derivative);
                    } else if (!rNode.IsFixed(instance_component)) {
                        rNode.FastGetSolutionStepValue(instance_component) =
                            previous_variable + this->GetDeltaTime() * previous_first_time_derivative +
                            0.5 * this->GetDeltaTime() * this->GetDeltaTime() * previous_second_time_derivative;
                    }
                }
            }
        });

        // Update (Angular) Acceleration, (Angular) Velocity and DtPressure
        this->UpdateVariablesDerivatives(rModelPart);

        KRATOS_CATCH("")
    }

protected:
    inline void UpdateVariablesDerivatives(ModelPart& rModelPart) override
    {
        KRATOS_TRY

        block_for_each(rModelPart.Nodes(), [this](Node& rNode) {
            // static here means no velocities and accelerations for the displacement/rotation D.O.F.
            for (const auto& r_first_order_scalar_variable : this->GetFirstOrderScalarVariables()) {
                this->UpdateScalarTimeDerivative(rNode, r_first_order_scalar_variable.instance,
                                                 r_first_order_scalar_variable.first_time_derivative);
            }
        });

        KRATOS_CATCH("")
    }

}; // Class NewmarkQuasistaticUPwScheme

} // namespace Kratos