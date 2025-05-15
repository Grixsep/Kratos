# TODO: this file does not seem to play any role ??

try:
    # import Kratos
    import KratosMultiphysics as KM

    # Import Kratos "wrapper" for unittests
    import KratosMultiphysics.KratosUnittest as KratosUnittest

    # Import the tests o test_classes to create the suits:
    import SmallTests
    import NightTests

except ImportError:
    print("Import error in test_LaserDrillingApplication.py")

KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)


def AssembleTestSuites():
    """Populates the test suites to run.

    Populates the test suites to run. At least, it should populate the suites:
    "small", "nightly" and "all"

    Return
    ------

    suites: A dictionary of suites
        The set of suites with its test_cases added.
    """
    # Suites to run
    suites = KratosUnittest.KratosSuites

    # SMALL TESTS
    small_suite = SmallTests.SetTestSuite(suites)

    # NIGHTLY TESTS
    night_suite = NightTests.SetTestSuite(suites)

    # include small suite in night suite
    night_suite.addTests(small_suite)

    # ALL TESTS
    all_suite = suites["all"]

    all_suite.addTests(night_suite)

    return suites


if __name__ == "__main__":
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.runTests(AssembleTestSuites())
