//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "droplet_dynamics_application.h"
#include "droplet_dynamics_application_variables.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"


namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosDropletDynamicsApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosDropletDynamicsApplication,
        KratosDropletDynamicsApplication::Pointer,
        KratosApplication>(m, "KratosDropletDynamicsApplication")
        .def(py::init<>())
        ;

    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,EXT_INT_FORCE)

    // Smoothed surface to calculate DISTANCE_GRADIENT
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DISTANCE_AUX);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,DISTANCE_AUX2);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,DISTANCE_GRADIENT_AUX);

    // Parallel levelset distance calculator needs an AREA_VARIABLE_AUX
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m,AREA_VARIABLE_AUX);

    // A variable to check if node is on cut element (maybe in a layer farther for future!)
    //KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, IS_NEAR_CUT)

    // Contact line calculation
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,NORMAL_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,TANGENT_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,CONTACT_VECTOR);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_ANGLE);
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m,CONTACT_VECTOR_MICRO);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_ANGLE_MICRO);
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONTACT_VELOCITY);

    // Enriched pressure is an array of NumNodes components defined for elements. Access it using Element.GetValue()
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ENRICHED_PRESSURE_1)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ENRICHED_PRESSURE_2)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ENRICHED_PRESSURE_3)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ENRICHED_PRESSURE_4)

    // Last known velocity and pressure to recalculate the last increment
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VELOCITY_STAR)
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PRESSURE_STAR)

    // Pressure gradient to calculate its jump over interface
    // KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, PRESSURE_GRADIENT_AUX)

    // Level-set convective velocity
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, CONVECTIVE_VELOCITY)
}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
