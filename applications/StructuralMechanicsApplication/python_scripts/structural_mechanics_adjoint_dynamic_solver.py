from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import MechanicalSolver

def CreateSolver(model, custom_settings):
    return StructuralMechanicsAdjointDynamicSolver(model, custom_settings)

class StructuralMechanicsAdjointDynamicSolver(MechanicalSolver):

    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[AdjointMechanicalSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "response_function_settings" : {},
            "sensitivity_settings" : {}
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def AddVariables(self):
        super().AddVariables()
        self._add_dynamic_variables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SHAPE_SENSITIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_VECTOR_2)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ADJOINT_VECTOR_3)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.AUX_ADJOINT_VECTOR_1)
        if self.settings["rotation_dofs"].GetBool():
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ADJOINT_ROTATION)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ANGULAR_ADJOINT_VECTOR_2)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ANGULAR_ADJOINT_VECTOR_3)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.ANGULAR_AUX_ADJOINT_VECTOR_1)
        # TODO evaluate if these variables should be historical
        KratosMultiphysics.Logger.PrintInfo("::[AdjointMechanicalSolver]:: ", "Variables ADDED")

    def PrepareModelPart(self):
        if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]!= 3):
            raise Exception("there are currently only 3D adjoint elements available")
        super().PrepareModelPart()
        # TODO Why does replacement need to happen after reading materials?

        process_info = self.main_model_part.ProcessInfo
        if (process_info.Has(StructuralMechanicsApplication.IS_ADJOINT) and
            process_info.GetValue(StructuralMechanicsApplication.IS_ADJOINT)):
            raise RuntimeError("Modelpart '{}' is already adjoint modelpart!".format(self.main_model_part.Name))

        # defines how the primal elements should be replaced with their adjoint counterparts
        replacement_settings = KratosMultiphysics.Parameters("""
            {
                "element_name_table" :
                {
                    "ShellThinElement3D3N"               : "AdjointFiniteDifferencingShellThinElement3D3N",
                    "CrLinearBeamElement3D2N"            : "AdjointFiniteDifferenceCrBeamElementLinear3D2N",
                    "TrussLinearElement3D2N"             : "AdjointFiniteDifferenceTrussLinearElement3D2N",
                    "TrussElement3D2N"                   : "AdjointFiniteDifferenceTrussElement3D2N",
                    "TotalLagrangianElement2D3N"         : "TotalLagrangianAdjointElement2D3N",
                    "TotalLagrangianElement2D4N"         : "TotalLagrangianAdjointElement2D4N",
                    "TotalLagrangianElement2D6N"         : "TotalLagrangianAdjointElement2D6N",
                    "TotalLagrangianElement3D4N"         : "TotalLagrangianAdjointElement3D4N",
                    "TotalLagrangianElement3D8N"         : "TotalLagrangianAdjointElement3D8N",
                    "SmallDisplacementElement3D4N"       : "AdjointFiniteDifferencingSmallDisplacementElement3D4N",
                    "SmallDisplacementElement3D6N"       : "AdjointFiniteDifferencingSmallDisplacementElement3D6N",
                    "SmallDisplacementElement3D8N"       : "AdjointFiniteDifferencingSmallDisplacementElement3D8N",
                    "SpringDamperElement3D2N"            : "AdjointFiniteDifferenceSpringDamperElement3D2N",
                    "NodalConcentratedDampedElement3D1N" : "AdjointFiniteDifferenceNodalConcentratedDampedElement3D1N"
                },
                "condition_name_table" :
                {
                    "PointLoadCondition2D1N"         : "AdjointSemiAnalyticPointLoadCondition2D1N",
                    "PointLoadCondition3D1N"         : "AdjointSemiAnalyticPointLoadCondition3D1N",
                    "SurfaceLoadCondition3D3N"       : "AdjointSemiAnalyticSurfaceLoadCondition3D3N",
                    "SurfaceLoadCondition3D4N"       : "AdjointSemiAnalyticSurfaceLoadCondition3D4N",
                    "SmallDisplacementSurfaceLoadCondition3D3N" : "AdjointSemiAnalyticSmallDisplacementSurfaceLoadCondition3D3N",
                    "SmallDisplacementSurfaceLoadCondition3D4N" : "AdjointSemiAnalyticSmallDisplacementSurfaceLoadCondition3D4N",
                    "LineLoadCondition3D2N"                     : "AdjointSemiAnalyticLineLoadCondition3D2N",
                    "SmallDisplacementLineLoadCondition3D2N"    : "AdjointSemiAnalyticSmallDisplacementLineLoadCondition3D2N"
                },
                "ignore_conditions" : [
                    "SurfaceCondition3D3N",
                    "SurfaceCondition3D4N",
                    "PointCondition3D1N"
                ]
            }
        """) # TODO remove "Condition3D" after issue#4439 is resolved

        StructuralMechanicsApplication.ReplaceMultipleElementsAndConditionsProcess(self.main_model_part, replacement_settings).Execute()
        process_info.SetValue(StructuralMechanicsApplication.IS_ADJOINT, True)

        KratosMultiphysics.Logger.PrintInfo("::[AdjointMechanicalSolver]:: ", "ModelPart prepared for Solver.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_DISPLACEMENT_Z, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ADJOINT_VECTOR_2_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ADJOINT_VECTOR_2_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ADJOINT_VECTOR_2_Z, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ADJOINT_VECTOR_3_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ADJOINT_VECTOR_3_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ADJOINT_VECTOR_3_Z, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.AUX_ADJOINT_VECTOR_1_X, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.AUX_ADJOINT_VECTOR_1_Y, self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.AUX_ADJOINT_VECTOR_1_Z, self.main_model_part)
        if self.settings["rotation_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_X, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Y, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ADJOINT_ROTATION_Z, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ANGULAR_ADJOINT_VECTOR_2_X, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ANGULAR_ADJOINT_VECTOR_2_Y, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ANGULAR_ADJOINT_VECTOR_2_Z, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ANGULAR_ADJOINT_VECTOR_3_X, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ANGULAR_ADJOINT_VECTOR_3_Y, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ANGULAR_ADJOINT_VECTOR_3_Z, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ANGULAR_AUX_ADJOINT_VECTOR_1_X, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ANGULAR_AUX_ADJOINT_VECTOR_1_Y, self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(StructuralMechanicsApplication.ANGULAR_AUX_ADJOINT_VECTOR_1_Z, self.main_model_part)            
        KratosMultiphysics.Logger.PrintInfo("::[AdjointMechanicalSolver]:: ", "DOF's ADDED.")

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        response_type = self.settings["response_function_settings"]["response_type"].GetString()
        if response_type == "adjoint_least_squares":
            # import measurement data
            if not self.settings["response_function_settings"].Has("measurement_data_filename"):
                raise Exception("AdjointMechanicalSolver: file with measurement data needed.")
            
            with open(self.settings["response_function_settings"]["measurement_data_filename"].GetString(), 'r') as parameter_file:
                measurement_data = KratosMultiphysics.Parameters(parameter_file.read())

            self.response_function = StructuralMechanicsApplication.AdjointLeastSquaresResponseFunction(self.main_model_part, self.settings["response_function_settings"], measurement_data)
        elif response_type == "adjoint_nodal_root_mean_square":
            self.response_function = StructuralMechanicsApplication.AdjointNodalRootMeanSquareResponseFunction(self.main_model_part, self.settings["response_function_settings"])
        elif response_type == "adjoint_damping_energy_dissipation":
            self.response_function = StructuralMechanicsApplication.AdjointDampingEnergyDissipationResponseFunction(self.main_model_part, self.settings["response_function_settings"])
        elif response_type == "adjoint_test":
            self.response_function = StructuralMechanicsApplication.AdjointTestResponseFunction(self.main_model_part, self.settings["response_function_settings"])
        else:
            raise Exception("invalid response_type: " + response_type)

        self.sensitivity_builder = KratosMultiphysics.SensitivityBuilder(self.settings["sensitivity_settings"], self.main_model_part, self.response_function)
        self.sensitivity_builder.Initialize()

        super().Initialize()
        self.response_function.Initialize()

        KratosMultiphysics.Logger.PrintInfo("::[AdjointMechanicalSolver]:: ", "Finished initialization.")

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        self.response_function.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.response_function.FinalizeSolutionStep()
        self.sensitivity_builder.UpdateSensitivities()

    def _create_mechanical_solution_strategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            if self.settings["compute_reactions"].GetBool():
                raise Exception("\"compute_reactions\" is not possible for adjoint models parts")
            if self.settings["move_mesh_flag"].GetBool():
                raise Exception("\"move_mesh_flag\" is not allowed for adjoint models parts")
            mechanical_solution_strategy = self._create_linear_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available for adjoints!\n"
            err_msg += "Available options are: \"linear\""
            raise Exception(err_msg)
        return mechanical_solution_strategy

    def _create_solution_scheme(self):
        scheme_type = self.settings["scheme_settings"]["scheme_type"].GetString()
        if scheme_type == "bossak":
            scheme = KratosMultiphysics.ResidualBasedAdjointBossakScheme(
                self.settings["scheme_settings"],
                self.response_function)
        else:
            raise Exception("Invalid scheme_type: " + scheme_type)
        return scheme
