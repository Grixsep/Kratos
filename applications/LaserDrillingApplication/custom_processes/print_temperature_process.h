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

#ifndef KRATOS_PRINT_TEMPERATURE_PROCESS_H
#define KRATOS_PRINT_TEMPERATURE_PROCESS_H

#include <pybind11/pybind11.h>

// System includes
#include <string>
#include <iostream>
#include <list>
// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"

// Application includes


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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

/// Utility to modify the distances of an embedded object in order to avoid bad intersections
/// Besides, it also deactivate the full negative distance elements
class KRATOS_API(LASER_DRILLING_APPLICATION) PrintTemperatureProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DistanceModificationProcess
    KRATOS_CLASS_POINTER_DEFINITION(PrintTemperatureProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Kratos parameters.
    PrintTemperatureProcess(
        ModelPart& rModelPart,
        Parameters& rParameters);

    /// Constructor with Kratos model
    PrintTemperatureProcess(
        Model& rModel,
        Parameters& rParameters);

    /// Destructor.
    ~PrintTemperatureProcess() override {}

    ///@}
    ///@name Operators
    ///@{

    void MakeMeasurements(pybind11::list &id_nodes,
                          pybind11::list &x_coord,
                          pybind11::list &y_coord,
                          pybind11::list &temperature);

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "PrintTemperatureProcess" ;
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {rOStream << "PrintTemperatureProcess";}

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {}


    ///@}
    ///@name Friends
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart&                                       mrModelPart;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    void CheckDefaultsAndProcessSettings(Parameters &rParameters);

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Default constructor.
    PrintTemperatureProcess() = delete;

    /// Assignment operator.
    PrintTemperatureProcess& operator=(PrintTemperatureProcess const& rOther) = delete;

    /// Copy constructor.
    PrintTemperatureProcess(PrintTemperatureProcess const& rOther) = delete;

    ///@}

}; // Class PrintTemperatureProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

///@} addtogroup block

};  // namespace Kratos.

#endif // KRATOS_PRINT_TEMPERATURE_PROCESS_H