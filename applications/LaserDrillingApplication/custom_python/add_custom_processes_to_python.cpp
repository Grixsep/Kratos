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

// External includes

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "includes/kratos_parameters.h"

#include "custom_processes/element_deactivation_process.hpp"
#include "custom_processes/print_temperature_process.h"
#include <pybind11/pybind11.h>

namespace Kratos
{

namespace Python
{

void AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<ElementDeactivationProcess, ElementDeactivationProcess::Pointer, Process>
    (m, "ElementDeactivationProcess", py::module_local())
    .def(py::init<ModelPart&, Parameters>());

    py::class_<PrintTemperatureProcess, PrintTemperatureProcess::Pointer, Process>
    (m, "PrintTemperatureProcess", py::module_local())
    .def(py::init<ModelPart&, Parameters&>())
    .def("MakeMeasurements", &PrintTemperatureProcess::MakeMeasurements);
}

}  // namespace Python.

} // Namespace Kratos
