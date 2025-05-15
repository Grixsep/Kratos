import os

# The aim of this test is to check that the temperature is positive at the end of each solution step

# Importing the Kratos Library
try:
    import KratosMultiphysics
except ImportError:
    print("Failed import of KratosMultiphysics")

# Import KratosUnittest
try:
    import KratosMultiphysics.KratosUnittest as KratosUnittest
except ImportError:
    print("Failed to import: import KratosMultiphysics.KratosUnittest as KratosUnittest")

try:
    from KratosMultiphysics.LaserDrillingApplication.laserdrilling_analysis import LaserDrillingAnalysis
except ImportError:
    print(
        "Failed to import: from KratosMultiphysics.LaserDrillingApplication.laser_drilling_analysis import LaserDrillingAnalysis"
    )

# This utility will control the execution scope

debug_mode = False


class controlledExecutionScope:
    def __init__(self, scope):
        self.currentPath = os.getcwd()
        self.scope = scope

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)


class PositiveTemperatureTestFactory(KratosUnittest.TestCase):
    def setUp(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            # Setting parameters

            with open(self.file_parameters, "r") as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())

            # Create Model
            model = KratosMultiphysics.Model()

            self.test = PositiveTemperatureAnalysis(model, parameters)

    def test_execution(self):
        with controlledExecutionScope(os.path.dirname(os.path.realpath(__file__))):
            self.test.Run()

    def tearDown(self):
        pass


class PositiveTemperatureAnalysis(LaserDrillingAnalysis):
    def __init__(self, model, parameters):
        super().__init__(model, parameters)
        if not debug_mode:
            parameters["output_processes"] = KratosMultiphysics.Parameters("""{}""")
        # self._GetLaserDrillingAnalysis().mdpas_folder_path = os.path.join(
        #     self._GetLaserDrillingAnalysis().main_path, "porosity_tests/porosity_conservation/"
        # )

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.CheckPositiveTemperature()

    def CheckPositiveTemperature(self):
        model_part = self.model.GetModelPart("ThermalModelPart")
        for node in model_part.Nodes:
            assert node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE) > 0
