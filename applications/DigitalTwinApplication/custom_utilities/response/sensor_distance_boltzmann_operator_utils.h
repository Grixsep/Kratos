//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "expression/container_expression.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class SensorDistanceBoltzmannOperatorResponseUtils
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensorDistanceBoltzmannOperatorResponseUtils);

    ///@}
    ///@name Life cycle
    ///@{

    SensorDistanceBoltzmannOperatorResponseUtils(
        ModelPart& rSensorModelPart,
        const double Beta);

    ///@}
    ///@name Public operations
    ///@{

    void Initialize();

    double CalculateValue();

    ContainerExpression<ModelPart::NodesContainerType> CalculateGradient() const;

    ///@}

private:
    ///@name Private member variables
    ///@{

    ModelPart* mpSensorModelPart;

    double mBeta;

    double mScaling;

    double mNumerator;

    double mDenominator;

    Matrix mDistances;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/