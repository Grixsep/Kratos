//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "utilities/openmp_utils.h"
#include "utilities/parallel_utilities.h"
#include "processes/find_nodal_h_process.h"

// Application includes
#include "print_temperature_process.h"
#include "laserdrilling_application_variables.h"

namespace Kratos
{


PrintTemperatureProcess::PrintTemperatureProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings

}

PrintTemperatureProcess::PrintTemperatureProcess(
    Model &rModel,
    Parameters &rParameters)
    : Process(),
      mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    // Check default settings

}

void PrintTemperatureProcess::CheckDefaultsAndProcessSettings(Parameters &rParameters)
{
    Parameters default_parameters( R"(
    {
    }  )" );

    rParameters.ValidateAndAssignDefaults(default_parameters);
}

void PrintTemperatureProcess::MakeMeasurements(pybind11::list &id_nodes, pybind11::list &x_coord, pybind11::list &y_coord, pybind11::list &temperature)
{
    KRATOS_TRY;

    const unsigned int n_elems = mrModelPart.NumberOfElements();
    std::vector<int> vector_of_id_nodes;

    for (unsigned int i_elem = 0; i_elem < n_elems; ++i_elem) {
        auto ind = mrModelPart.ElementsBegin() + i_elem;
        if (ind->Is(ACTIVE)){
            const auto NumNodes = ind->GetGeometry().PointsNumber();
            for (unsigned int j = 0; j < NumNodes; ++j){
                unsigned int id_node = ind->GetGeometry()[j].Id();
                if (std::find(vector_of_id_nodes.begin(), vector_of_id_nodes.end(), id_node) == vector_of_id_nodes.end()) {
                    vector_of_id_nodes.push_back(id_node);
                    id_nodes.append(id_node);
                    x_coord.append(ind->GetGeometry()[j].X());
                    y_coord.append(ind->GetGeometry()[j].Y());
                    temperature.append(ind->GetGeometry()[j].GetSolutionStepValue(TEMPERATURE));
                }
            }
        }
    }

    KRATOS_CATCH("");
}

/* Private functions ****************************************************/

};  // namespace Kratos.