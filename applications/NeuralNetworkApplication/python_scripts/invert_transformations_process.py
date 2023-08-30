import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess
from KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities import ImportDataFromFile
import numpy as np


def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a parameters object, encapsulating a json string")
    return InvertTransformationsProcess(settings["parameters"])

class InvertTransformationsProcess(NeuralNetworkProcess):

    def __init__(self,parameters):

        super().__init__()


        default_settings = KM.Parameters("""{
            "echo"                : 0,
            "input_file"          : "",
            "predictions_file"    : "",
            "input_file_name"     : "",
            "output_file_name"    : "",
            "write_output"        : true,
            "to_output"           : false,
            "save_format"         : "ascii"            
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.echo = parameters["echo"].GetInt()
        self.input_file = parameters["input_file"].GetString()
        self.predictions_file = parameters["predictions_file"].GetString()
        self.input_file_name = parameters["input_file_name"].GetString()
        self.output_file_name = parameters["output_file_name"].GetString()
        self.to_output = parameters["to_output"].GetBool()
        self.write_output = parameters["write_output"].GetBool()
        self.save_format = parameters["save_format"].GetString()
        if  not (self.save_format == "ascii") and  not (self.save_format == "npy"):
            raise Exception("Saving format not supported. Supported formats are ascii and npy.")

    def TransformPredictions(self, processes, data_in = None, data_out = None):
        # Load the data
        if data_in is None:
            self.data_in = ImportDataFromFile(self.input_file, "InputData")
        else:
            self.data_in = data_in
        if data_out is None:
            self.data_out = ImportDataFromFile(self.predictions_file, "OutputData")
        else:
            self.data_out = data_out

        # Invert the transformations
        for process in reversed(processes):
            try:
                inversion_settings = process.Invert(self.data_in, self.data_out)
                if not inversion_settings is None:
                    self.data_in = inversion_settings[0]
                    self.data_out = inversion_settings[1]
            except AttributeError:
                pass

        # Save the data
        if self.write_output:
            if self.save_format == "ascii":
                extension  = ".csv"
                self.input_file_name = self.input_file_name + extension
                self.output_file_name = self.output_file_name + extension
                np.savetxt(self.input_file_name, self.data_in.ExportAsArray())
                np.savetxt(self.output_file_name,self.data_out.ExportAsArray())

            elif self.save_format == "npy":
                extension = ".npy"
                self.input_file_name = self.input_file_name + extension
                self.output_file_name = self.output_file_name + extension
                np.save(self.input_file_name, self.data_in.ExportAsArray())
                np.save(self.output_file_name,self.data_out.ExportAsArray())
            else:
                raise Exception("Output save format not supported. Supported formats are ascii and npy.")

        # Return
        if self.to_output:
            return self.data_out
