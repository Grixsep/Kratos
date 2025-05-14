# Definition of the classes for the SMALL TESTS

# Import Kratos and necessary applications
import KratosMultiphysics
import KratosMultiphysics.LaserDrillingApplication
from KratosMultiphysics import Logger

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as UnitTest

# Importing test factories if possible
try:
    import TemperatureTestFactory as TemperatureTF

    temperature_imports_available = True
except ImportError:
    temperature_imports_available = False


class temperature_ranges_test(TemperatureTF.TemperatureRangesTestFactory):
    file_name = "porosity_tests/porosity_conservation/Test_porosityFluid"
    file_parameters = "porosity_tests/porosity_conservation/ProjectParameters.json"


available_tests = []
available_tests += [test_class for test_class in PorosityTF.PorosityConservationTestFactory.__subclasses__()]


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
