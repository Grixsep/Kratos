//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//                  Suneth Warnakulasuriya
//

// System includes

// Project includes
#include "utilities/data_type_traits.h"

// Application includes

// Include base h
#include "custom_io/hdf5_data_value_container_io.h"

namespace Kratos
{
namespace DataValueContainerIOHelperUtilities
{

template<class TDataType>
bool ReadComponentData(
    HDF5::File& rFile,
    const std::string& rVariableName,
    const std::string& rPrefix,
    DataValueContainer& rData)
{
    KRATOS_TRY

    if (KratosComponents<Variable<TDataType>>::Has(rVariableName)) {
        const auto& r_variable = KratosComponents<Variable<TDataType>>::Get(rVariableName);
        TDataType value{};
        rFile.ReadAttribute(rPrefix + "/DataValues", rVariableName, value);
        rData[r_variable] = value;
        return true;
    } else {
        return false;
    }

    KRATOS_CATCH("");
}

template<class TDataType>
bool WriteComponentData(
    HDF5::File& rFile,
    const std::string& rVariableName,
    const std::string& rPrefix,
    const DataValueContainer& rData)
{
    KRATOS_TRY

    if (KratosComponents<Variable<TDataType>>::Has(rVariableName)) {
        const auto& r_variable = KratosComponents<Variable<TDataType>>::Get(rVariableName);
        const auto& r_value = rData[r_variable];
        rFile.WriteAttribute(rPrefix + "/DataValues", rVariableName, r_value);
        return true;
    } else {
        return false;
    }

    KRATOS_CATCH("");
}

template<class... TDataTypes>
void Read(
    HDF5::File& rFile,
    const std::string& rVariableName,
    const std::string& rPrefix,
    DataValueContainer& rData)
{
    KRATOS_TRY

    const bool is_read = (... || ReadComponentData<TDataTypes>(rFile, rVariableName, rPrefix, rData));

    // KRATOS_ERROR_IF_NOT(is_read) << "The variable \"" << rVariableName << "\" not found in registered variables list.";

    KRATOS_CATCH("");
}

template<class... TDataTypes>
void Write(
    HDF5::File& rFile,
    const std::string& rVariableName,
    const std::string& rPrefix,
    const DataValueContainer& rData)
{
    KRATOS_TRY

    const bool is_written = (... || WriteComponentData<TDataTypes>(rFile, rVariableName, rPrefix, rData));

    // KRATOS_ERROR_IF_NOT(is_written) << "The variable \"" << rVariableName << "\" not found in registered variables list.";

    KRATOS_CATCH("");
}

} // namespace DataValueContainerIOHelperUtilities

namespace HDF5
{
namespace Internals
{

void ReadDataValueContainer(
    File& rFile,
    const std::string& rPrefix,
    DataValueContainer& rData)
{
    KRATOS_TRY;

    const auto& attr_names = rFile.GetAttributeNames(rPrefix + "/DataValues");

    for (const auto& r_name : attr_names) {
        DataValueContainerIOHelperUtilities::Read<
            int,
            double,
            std::string,
            array_1d<double, 3>,
            array_1d<double, 4>,
            array_1d<double, 6>,
            array_1d<double, 9>,
            Kratos::Vector,
            Kratos::Matrix>(rFile, r_name, rPrefix, rData);
    }

    KRATOS_CATCH("Path: \"" + rPrefix + "/DataValues\".");
}

void WriteDataValueContainer(
    File& rFile,
    const std::string& rPrefix,
    const DataValueContainer& rData)
{
    KRATOS_TRY;

    rFile.AddPath(rPrefix + "/DataValues");

    for (auto it = rData.begin(); it != rData.end(); ++it) {
        DataValueContainerIOHelperUtilities::Write<
            int,
            double,
            std::string,
            array_1d<double, 3>,
            array_1d<double, 4>,
            array_1d<double, 6>,
            array_1d<double, 9>,
            Kratos::Vector,
            Kratos::Matrix>(rFile, it->first->Name(), rPrefix, rData);
    }

    KRATOS_CATCH("Path: \"" + rPrefix + "/DataValues\".");
}

} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.