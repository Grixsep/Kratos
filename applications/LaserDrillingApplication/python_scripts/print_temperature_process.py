import os

import numpy as np
import h5py

import KratosMultiphysics as Kratos
import KratosMultiphysics.LaserDrillingApplication as KratosLaserDrilling


def Factory(settings, Model):
    if type(settings) != Kratos.Parameters:  # TODO: Ruff E721
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return PrintTemperatureProcess(Model, settings["Parameters"])

class PrintTemperatureProcess(Kratos.Process):

    def __init__(self, Model, settings):

        Kratos.Process.__init__(self)

        default_parameters = Kratos.Parameters( """
        {
            "model_part_name"                : "CHOOSE_FLUID_MODELPART_NAME",
            "file_name"                      : "temperature_db.hdf5",
            "step_jump_between_outputs"      : 1
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters)

        self.problem_path = os.getcwd()
        self.file_path = os.path.join(str(self.problem_path),settings["file_name"].GetString())
        self.step_jump_between_outputs = settings["step_jump_between_outputs"].GetInt()
        self.main_model_part = Model[settings["model_part_name"].GetString()]
        self.print_number = 0
        self.jump_between_outputs_counter = 0

    def ExecuteInitialize(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        self.jump_between_outputs_counter += 1
        if self.jump_between_outputs_counter >= self.step_jump_between_outputs:
            self.jump_between_pulses_counter = 0
            self.x_coord = []
            self.y_coord = []
            self.id_node = []
            self.temperature = []
            self.print_number += 1
            self.time = self.main_model_part.ProcessInfo[Kratos.TIME]

            for elem in self.main_model_part.Elements:
                if elem.Is(Kratos.ACTIVE):
                    for node in elem.GetNodes():
                        if not node.Id in self.id_node:
                            self.id_node.append(node.Id)
                            self.x_coord.append(node.X)
                            self.y_coord.append(node.Y)
                            self.temperature.append(node.GetSolutionStepValue(Kratos.TEMPERATURE))

            with h5py.File(self.file_path, 'a') as f:
                    self.WriteDataToFile(file_or_group = f,
                                names = ['Id', 'x_coord', 'y_coord', 'temperature'],
                                data = [self.id_node, self.x_coord, self.y_coord, self.temperature])

    def WriteDataToFile(self, file_or_group, names, data):
        if str(self.print_number) in file_or_group:
            file_or_group['/'].__delitem__(str(self.print_number))
        self.sub_group = file_or_group.create_group(str(self.print_number))
        self.sub_group.attrs['time'] = str(self.time)

        column_name = np.dtype({'names': names, 'formats':[(float)]*len(names)})
        data_array = np.rec.fromarrays(data, dtype = column_name)

        self.sub_group.create_dataset(name = str(self.print_number), data = data_array)