# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.SwimmingDEMApplication as KratosSDEM

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplySkrzypaczPorositySolutionAndBodyForceProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class ApplySkrzypaczPorositySolutionAndBodyForceProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings):
        """The default constructor of the class.

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the fluid model part.
        settings -- Kratos parameters containing process settings.
        """
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
            {
                "model_part_name"          : "please_specify_model_part_name",
                "variable_name"            : "BODY_FORCE",
                "benchmark_name"           : "custom_body_force.vortex",
                "benchmark_parameters"     : {},
                "compute_nodal_error"      : true,
                "print_convergence_output" : false,
                "output_parameters"        : {}
            }
            """
            )

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[self.settings["model_part_name"].GetString()]
        self.variable = KratosMultiphysics.KratosGlobals.GetVariable(self.settings["variable_name"].GetString())

        self.ApplySkrzypaczPorositySolutionAndBodyForceProcess = KratosSDEM.SkrzypaczPorositySolutionAndBodyForceProcess(self.model_part, settings)


    def ExecuteBeforeSolutionLoop(self):
        self.ApplySkrzypaczPorositySolutionAndBodyForceProcess.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.ApplySkrzypaczPorositySolutionAndBodyForceProcess.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        pass