# Definition of the classes for the SMALL TESTS

# Import Kratos and necessary applications
from KratosMultiphysics import Logger


# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as UnitTest

# Importing test factories if possible
try:
    import TemperatureTestFactory as TemperatureTF

    temperature_imports_available = True
except ImportError:
    Logger.PrintWarning("Failed to import TemperatureTestFactory")
    temperature_imports_available = False


class positive_temperature_test(TemperatureTF.PositiveTemperatureTestFactory):
    file_name = "temperature_tests/positive_temperature/Test_positive_temperature"  # MDPA filename
    file_parameters = (
        "temperature_tests/positive_temperature/parameters/ProjectParameters.json"  # ProjectParameters filename
    )


available_tests = []
available_tests += [test_class for test_class in TemperatureTF.PositiveTemperatureTestFactory.__subclasses__()]


def SetTestSuite(suites):
    small_suite = suites["small"]
    small_suite.addTests(UnitTest.TestLoader().loadTestsFromTestCases(available_tests))

    return small_suite


def AssembleTestSuites():
    suites = UnitTest.KratosSuites
    small_suite = SetTestSuite(suites)
    suites["all"].addTests(small_suite)
    return suites


if __name__ == "__main__":
    debug_mode = False
    if debug_mode:
        severity = Logger.Severity.DETAIL
    else:
        severity = Logger.Severity.WARNING
    Logger.GetDefaultOutput().SetSeverity(severity)
    UnitTest.runTests(AssembleTestSuites())
    UnitTest.main()
