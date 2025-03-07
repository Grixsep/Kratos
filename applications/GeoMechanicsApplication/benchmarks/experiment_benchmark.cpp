//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Richard Faasse
//

#include <benchmark/benchmark.h>

#include "custom_elements/plane_strain_stress_state.h"
#include "custom_elements/three_dimensional_stress_state.h"
#include "custom_elements/transient_Pw_element.hpp"
#include <boost/numeric/ublas/assignment.hpp>
#include "includes/model_part.h"
#include "containers/model.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "includes/cfd_variables.h"
#include "geo_mechanics_application_variables.h"
namespace
{
using namespace Kratos;

PointerVector<Node> CreateThreeNodes()
{
    PointerVector<Node> result;
    result.push_back(make_intrusive<Node>(1, 0.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(2, 1.0, 0.0, 0.0));
    result.push_back(make_intrusive<Node>(3, 1.0, 1.0, 0.0));
    return result;
}

PointerVector<Node> CreateThreeCoincidentNodes()
{
    PointerVector<Node> result;
    for (unsigned int id = 1; id <= 3; id++) {
        result.push_back(make_intrusive<Node>(id, 0.0, 0.0, 0.0));
    }
    return result;
}

template <unsigned int TNumNodes>
PointerVector<Node> CreateNodesOnModelPart(ModelPart& rModelPart)
{
    PointerVector<Node> result;
    result.push_back(rModelPart.CreateNewNode(1, 0.0, 0.0, 0.0));
    result.push_back(rModelPart.CreateNewNode(2, 1.0, 0.0, 0.0));
    result.push_back(rModelPart.CreateNewNode(3, 1.0, 1.0, 0.0));
    if constexpr (TNumNodes == 4) {
        result.push_back(rModelPart.CreateNewNode(4, 1.0, 1.0, 1.0));
    }
    return result;
}

ModelPart& CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(Model& rModel)
{
    auto& r_result = rModel.CreateModelPart("Main");
    r_result.AddNodalSolutionStepVariable(WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(DT_WATER_PRESSURE);
    r_result.AddNodalSolutionStepVariable(VOLUME_ACCELERATION);

    return r_result;
}

void RemoveThreeNodes(ModelPart& rModelPart)
{
    rModelPart.RemoveNodeFromAllLevels(1);
    rModelPart.RemoveNodeFromAllLevels(2);
    rModelPart.RemoveNodeFromAllLevels(3);
}

Element::IndexType NextElementNumber(const ModelPart& rModelPart)
{
    return rModelPart.NumberOfElements() + 1;
}

template <unsigned int TDim, unsigned int TNumNodes>
intrusive_ptr<TransientPwElement<TDim, TNumNodes>> CreateTransientPwElementWithPWDofs(ModelPart& rModelPart,
                                                                                      const Properties::Pointer& rProperties)
{
    intrusive_ptr<TransientPwElement<TDim, TNumNodes>> p_element;
    if constexpr (TDim == 2) {
        p_element = make_intrusive<TransientPwElement<TDim, TNumNodes>>(
            NextElementNumber(rModelPart),
            std::make_shared<Triangle2D3<Node>>(CreateNodesOnModelPart<TNumNodes>(rModelPart)),
            rProperties, std::make_unique<PlaneStrainStressState>());
    } else {
        p_element = make_intrusive<TransientPwElement<TDim, TNumNodes>>(
            NextElementNumber(rModelPart),
            std::make_shared<Tetrahedra3D4<Node>>(CreateNodesOnModelPart<TNumNodes>(rModelPart)),
            rProperties, std::make_unique<ThreeDimensionalStressState>());
    }
    for (auto& r_node : p_element->GetGeometry()) {
        r_node.AddDof(WATER_PRESSURE);
    }
    rModelPart.AddElement(p_element);
    return p_element;
}

intrusive_ptr<TransientPwElement<2, 3>> CreateTriangleTransientPwElementWithoutPWDofs(ModelPart& rModelPart,
                                                                                      const Properties::Pointer& rProperties)
{
    auto p_element = make_intrusive<TransientPwElement<2, 3>>(
        NextElementNumber(rModelPart),
        std::make_shared<Triangle2D3<Node>>(CreateNodesOnModelPart<3>(rModelPart)), rProperties,
        std::make_unique<PlaneStrainStressState>());

    rModelPart.AddElement(p_element);
    return p_element;
}

template <unsigned int TDim, unsigned int TNumNodes>
void SetBasicPropertiesAndVariables(intrusive_ptr<TransientPwElement<TDim, TNumNodes>> rElement)
{
    rElement->GetProperties().SetValue(DENSITY_WATER, 1.0E3);
    rElement->GetProperties().SetValue(DYNAMIC_VISCOSITY, 1.0E-2);
    rElement->GetProperties().SetValue(PERMEABILITY_XX, 1.0);
    rElement->GetProperties().SetValue(PERMEABILITY_YY, 1.0);
    rElement->GetProperties().SetValue(PERMEABILITY_XY, 1.0);
    if constexpr (TDim == 3) {
        rElement->GetProperties().SetValue(PERMEABILITY_ZZ, 1.0);
        rElement->GetProperties().SetValue(PERMEABILITY_YZ, 1.0);
        rElement->GetProperties().SetValue(PERMEABILITY_ZX, 1.0);
    }
    const auto gravity_acceleration = array_1d<double, 3>{0.0, -10.0, 0.0};
    for (auto& r_node : rElement->GetGeometry()) {
        r_node.FastGetSolutionStepValue(VOLUME_ACCELERATION) = gravity_acceleration;
        r_node.FastGetSolutionStepValue(WATER_PRESSURE)      = 0.0;
    }
}

} // namespace

namespace Kratos
{

void testBenchmark(benchmark::State& state)
{
    // Arrange
    Model model;
    auto& r_model_part = CreateModelPartWithWaterPressureVariableAndVolumeAcceleration(model);
    r_model_part.AddNodalSolutionStepVariable(HYDRAULIC_DISCHARGE);
    auto p_element = CreateTransientPwElementWithPWDofs<2, 3>(r_model_part, std::make_shared<Properties>());
    SetBasicPropertiesAndVariables(p_element);
    p_element->GetProperties().SetValue(BIOT_COEFFICIENT, 0.5);
    p_element->GetProperties().SetValue(BULK_MODULUS_FLUID, 1.0E6);
    p_element->GetProperties().SetValue(BULK_MODULUS_SOLID, 1.0E6);
    p_element->GetProperties().SetValue(POROSITY, 0.1);
    p_element->GetProperties().SetValue(IGNORE_UNDRAINED, false);
    const auto dummy_process_info = ProcessInfo{};
    p_element->Initialize(dummy_process_info);
    p_element->InitializeSolutionStep(dummy_process_info);

    for (auto _ : state) {
        Vector actual_right_hand_side;
        Matrix actual_left_hand_side;
        p_element->CalculateLocalSystem(actual_left_hand_side, actual_right_hand_side, dummy_process_info);
    }
}


BENCHMARK(testBenchmark);

}  // namespace Kratos

BENCHMARK_MAIN();

