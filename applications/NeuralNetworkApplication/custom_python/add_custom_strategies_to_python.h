//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Rishith Ellath Meethal (https://github.com/rishithellathmeethal)
//                   Daniel Andrés Arcones https://github.com/danielandresarcones
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define_python.h"

namespace Kratos::Python {

void  AddCustomStrategiesToPython(pybind11::module& m);

}  // namespace Kratos::Python.
