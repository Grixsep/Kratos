import KratosMultiphysics as KM
import KratosMultiphysics.NeuralNetworkApplication.input_dataclasses as InputDataclasses
from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData
import numpy as np

# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing NeuralNetworkApplication
if not CheckIfApplicationsAvailable("NeuralNetworkApplication"):
    raise ImportError("The NeuralNetworkApplication is not available!")
from KratosMultiphysics.NeuralNetworkApplication.neural_network_analysis import NeuralNetworkAnalysis

def Create(settings, model, solver_name):
    return NeuralNetworkWrapper(settings, model, solver_name)

class NeuralNetworkWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the NeuralNetworkApplication of Kratos"""
    def __init__(self, settings, model, solver_name):
        super(NeuralNetworkWrapper, self).__init__(settings, model, solver_name)

        input_file_name = self.settings["solver_wrapper_settings"]["input_file"].GetString()
        # TODO: The input and output can only be one variable in the current state.
        # TODO: Only nodal variables can be used in the current state.
        self.input_variable = self.settings["solver_wrapper_settings"]["input_variable"].GetString()
        self.output_variable = self.settings["solver_wrapper_settings"]["output_variable"].GetString()
        self.input_variable_name = self.settings["solver_wrapper_settings"]["input_variable_name"].GetString()
        self.output_variable_name = self.settings["solver_wrapper_settings"]["output_variable_name"].GetString()
        self.time_buffer = self.settings["solver_wrapper_settings"]["time_buffer"].GetInt()
        # self.submodel_part_name = self.settings["solver_wrapper_settings"]["submodel_part_name"].GetString()
        self.mp = self.model.CreateModelPart("NeuralNetwork")
        # self.mp.CreateSubModelPart(self.submodel_part_name)
        self.mp.AddNodalSolutionStepVariable(getattr(KM, self.input_variable_name))
        self.mp.AddNodalSolutionStepVariable(getattr(KM, self.output_variable_name))
        mdpa_file_name = self.settings["solver_wrapper_settings"]["mdpa_file_name"].GetString()

        self.mp.ProcessInfo[KM.DOMAIN_SIZE] = 2

        self.lookback = self.settings["solver_wrapper_settings"] ["lookback"].GetInt()
        try:
            self.record = settings["solver_wrapper_settings"]["record"].GetBool()
        except RuntimeError:
            self.record = False
        try:
            self.only_input = settings["solver_wrapper_settings"]["only_input"].GetBool()
        except RuntimeError:
            self.only_input = False
        try:
            self.timesteps_as_features = settings["solver_wrapper_settings"]["timesteps_as_features"].GetBool()
        except RuntimeError:
            self.timesteps_as_features = False
        try:
            self.feaures_as_timestpes = settings["solver_wrapper_settings"]["features_as_timesteps"].GetBool()
        except RuntimeError:
            self.feaures_as_timestpes = False
        try:
            self.soft_start_flag = settings["solver_wrapper_settings"]["soft_start"].GetBool()
        except RuntimeError:
            self.soft_start_flag = True
        try:
            self.dimension_in = settings["solver_wrapper_settings"]["dimension_input"].GetInt()
        except RuntimeError:
            self.dimension_in = 1
        try:
            self.reorder_partitions = settings["solver_wrapper_settings"]["reorder_partitions"].GetInt()
        except RuntimeError:
            self.reorder_partitions = 1
    
        KM.ModelPartIO(mdpa_file_name).ReadModelPart(self.mp)

        # self.mp.ProcessInfo[KM.DOMAIN_SIZE] = self.project_parameters["domain_size"].GetInt()

    def AdvanceInTime(self, current_time):
        return 0.0

    def SolveSolutionStep(self):
        with self.thread_manager:
            print("Receiving data into the neural network model")

            # Restarting preprocessed data if no FinalizeSolutionStep took place
            try:
                if self.converge_flag == False:
                    self.preprocessed_data_structure.UpdateData(self.preprocessed_previous.data)
                    if hasattr(self.preprocessed_previous,'lookback_data'):
                        if not self.only_input:
                            self.preprocessed_data_structure.UpdateLookbackAll(self.preprocessed_previous.lookback_data)
                    elif self.preprocessed_data_structure.data is None:
                        if not self.only_input:
                            self.preprocessed_data_structure.UpdateLookbackAll(np.zeros_like(self.preprocessed_data_structure.lookback_data))
                        self.soft_start_flag = True
                        if self.only_input:
                            if self.record:
                                self.preprocessed_data_structure.record_data = False
                            self.preprocessed_data_structure.lookback_state = False
                            delattr(self.preprocessed_data_structure, 'lookback_data')
                    else:
                        if self.only_input:
                            delattr(self.preprocessed_data_structure, 'lookback_data')
                            self.preprocessed_data_structure.lookback_state = False
                        else:
                            self.preprocessed_data_structure.lookback_state = False
                            if self.record:
                                self.preprocessed_data_structure.record_data = False
            except AttributeError:
                pass

            # Retrieve input from interface
            from_interface = self.GetInterfaceData(self.input_variable).GetData()
            from_interface = np.reshape(from_interface, (int(from_interface.size/self.dimension_in), self.dimension_in))
            self.input_data_structure.UpdateData(from_interface)

            # Initialize output from interface in first iteration
            if self.output_data_structure.data is None:
                self.output_data_structure.UpdateData(self.GetInterfaceData(self.output_variable).GetData())
            
            # Preprocess input and output 
            [self.input_data_structure, self.output_data_structure] = self._analysis_stage.Preprocessing(
                data_in = self.input_data_structure, data_out = self.output_data_structure)
            # Initialize preprocessed input to the network in first iteration
            if self.lookback>0 and not self.preprocessed_data_structure.lookback_state:
                if not self.only_input:
                    self.preprocessed_data_structure.CheckLookbackAndUpdate(self.output_data_structure.ExportAsArray())
                if self.record:
                    self.preprocessed_data_structure.CheckRecordAndUpdate(self.input_data_structure.ExportAsArray())
            # Update input to the structure with the new preprocessed input
            if self.record:
                if self.preprocessed_data_structure.data is None:
                    self.preprocessed_data_structure.record_data = False
                    self.preprocessed_data_structure.CheckRecordAndUpdate(self.input_data_structure.ExportAsArray())
                # This flag softens the initial lookback in the record (stabilizes the behaviour)
                if self.soft_start_flag:
                    for i in range(self.lookback-1):
                        if not self.only_input:
                            self.preprocessed_data_structure.CheckLookbackAndUpdate(self.output_data_structure.ExportAsArray())
                        self.preprocessed_data_structure.UpdateRecordLast(self.input_data_structure.ExportAsArray())
                    self.preprocessed_data_structure.UpdateRecordLast(self.input_data_structure.ExportAsArray())
                    self.soft_start_flag = False
                else:
                    self.preprocessed_data_structure.UpdateRecordLast(self.input_data_structure.ExportAsArray())
                self.record_data = self.preprocessed_data_structure.data
            else:
                self.preprocessed_data_structure.UpdateData(self.input_data_structure.ExportAsArray())
            # Set the reordering of the preprocessed data structure
            if self.reorder_partitions > 1:
                self.preprocessed_data_structure.reorder_partitions = self.reorder_partitions
                
            # Predict (and invert the transformations) from the new input and update it to the output
            self.output_data_structure.UpdateData(self._analysis_stage.Predict(data_structure_in = self.preprocessed_data_structure))
            # Export output to the interface if the buffer is full
            if self.current_time >= self.time_buffer:
                self.GetInterfaceData(self.output_variable).SetData(np.squeeze(self.output_data_structure.ExportAsArray()))
            print("Predicted data using neural network model.")

            
            

    def Initialize(self):
        self.input_data_structure = InputDataclasses.NeuralNetworkData()
        if self.lookback>0:
            self.preprocessed_data_structure = InputDataclasses.DataWithLookback(lookback_index=self.lookback)  
            self.preprocessed_previous = InputDataclasses.DataWithLookback(lookback_index=self.lookback)
        else:
            self.preprocessed_data_structure = InputDataclasses.NeuralNetworkData()
            self.preprocessed_previous = InputDataclasses.NeuralNetworkData()
        self.output_data_structure = InputDataclasses.NeuralNetworkData()
        self.current_time = 0
        
        self.data_dict = {data_name : CouplingInterfaceData(data_config, self.model, data_name, self.name) for (data_name, data_config) in self.settings["data"].items()}

    def Finalize(self):
        pass

    def InitializeSolutionStep(self):
        self.converge_flag = False

    def Predict(self):
        pass

    def FinalizeSolutionStep(self):
        if self.record:
            self.preprocessed_data_structure.data=self.record_data
        self.preprocessed_previous.UpdateData(self.preprocessed_data_structure.data)
        if not self.only_input:
            self.preprocessed_previous.UpdateLookbackAll(self.preprocessed_data_structure.lookback_data)
            self.preprocessed_previous.reorder_partitions = self.reorder_partitions
            self.reorder_partitions = 0 # Flag for only setting the reorder once
        self.converge_flag = True
        self.current_time += 1

    def OutputSolutionStep(self):
        pass

    def _CreateAnalysisStage(self):
        return NeuralNetworkAnalysis(self.project_parameters)
