//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


// System includes

// External includes
//

// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "laserdrilling_application.h"
#include "includes/variables.h"

namespace Kratos {



KratosLaserDrillingApplication::KratosLaserDrillingApplication()
    : KratosApplication("LaserDrillingApplication"),
      mLaserAxisymmetricEulerianConvectionDiffusion2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mLaserAxisymmetricEulerianConvectionDiffusion2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node >(Element::GeometryType::PointsArrayType(4))))
      {}

void KratosLaserDrillingApplication::Register() {
    KRATOS_INFO("") << "Initializing KratosLaserDrillingApplication... " << std::endl;

    // Registering variables
    KRATOS_REGISTER_VARIABLE(THERMAL_ENERGY)
    KRATOS_REGISTER_VARIABLE(THERMAL_ENERGY_PER_VOLUME)
    KRATOS_REGISTER_VARIABLE(THERMAL_DECOMPOSITION)
    KRATOS_REGISTER_VARIABLE(DECOMPOSITION_LAW_REFERENCE_TEMPERATURE)
    KRATOS_REGISTER_VARIABLE(DECOMPOSITION_LAW_CONSTANT_1)
    KRATOS_REGISTER_VARIABLE(DECOMPOSITION_LAW_CONSTANT_2)
    KRATOS_REGISTER_VARIABLE(THERMAL_ENERGY_PER_VOLUME_THRESHOLD)
    KRATOS_REGISTER_VARIABLE(ALPHA_THRESHOLD)
    KRATOS_REGISTER_VARIABLE(THERMAL_COUNTER)
    KRATOS_REGISTER_VARIABLE(DECOMPOSITION_LAW)
    KRATOS_REGISTER_VARIABLE(DECOMPOSED_ELEMENTAL_VOLUME)

    // Registering elements and conditions here
    KRATOS_REGISTER_ELEMENT("LaserAxisymmetricEulerianConvectionDiffusion2D3N", mLaserAxisymmetricEulerianConvectionDiffusion2D3N);
    KRATOS_REGISTER_ELEMENT("LaserAxisymmetricEulerianConvectionDiffusion2D4N", mLaserAxisymmetricEulerianConvectionDiffusion2D4N);

}

}  // namespace Kratos.
