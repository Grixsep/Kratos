//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Marco Antonio Zuñiga Perez
//

// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_strategies_to_python.h"

#include "processes/process.h"
#include "custom_utilities/solver_settings.h"

#include "spaces/ublas_space.h"

//schemes
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "custom_strategies/transonic_coefficients_update_scheme.h"

//linear solvers
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
namespace Python
{

void AddCustomStrategiesToPython(pybind11::module &m)
{
    namespace py = pybind11;
    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double>> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType > ResidualBasedIncrementalUpdateStaticSchemeType;

    //********************************************************************
    typedef TransonicCoefficientsUpdateScheme<SparseSpaceType, LocalSpaceType> TransonicCoefficientsUpdateSchemeType;
    py::class_<TransonicCoefficientsUpdateSchemeType, typename TransonicCoefficientsUpdateSchemeType::Pointer, ResidualBasedIncrementalUpdateStaticSchemeType>(m, "TransonicCoefficientsUpdateScheme")
        .def(py::init<Parameters >() )
        .def(py::init< >()
        );
}

} // namespace Python.

} // Namespace Kratos
