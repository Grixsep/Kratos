import KratosMultiphysics
from KratosMultiphysics.LaserDrillingApplication.laserdrilling_analysis import LaserDrillingAnalysis
import os

if __name__ == "__main__":
    from sys import argv

    # Define the directory containing the JSON files
    parameter_dir = "parameter_jsons"

    # Check for correct number of input arguments
    if len(argv) > 2:
        err_msg = "Too many input arguments!\n"
        err_msg += "Use this script in the following way:\n"
        err_msg += '- With default ProjectParameters directory ("parameter_jsons"):\n'
        err_msg += '    "python3 laserdrilling_analysis.py"\n'
        err_msg += "- With custom ProjectParameters file:\n"
        err_msg += '    "python3 laserdrilling_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2:  # Single ProjectParameters file is being passed
        project_parameters_file_name = argv[1]
        file_list = [project_parameters_file_name]
    else:  # Load all JSON files from parameter_jsons directory
        if not os.path.exists(parameter_dir):
            raise Exception(f"Directory '{parameter_dir}' does not exist!")
        file_list = [
            os.path.join(parameter_dir, f)
            for f in os.listdir(parameter_dir)
            if f.endswith(".json") and os.path.isfile(os.path.join(parameter_dir, f))
        ]

    if not file_list:
        raise Exception(f"No JSON files found in directory '{parameter_dir}'!")

    # Run simulation for each JSON file
    for project_parameters_file_name in file_list:
        print(f"------------ Running simulation with parameters from: {project_parameters_file_name} -------------")

        with open(project_parameters_file_name, "r") as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        model = KratosMultiphysics.Model()
        simulation = LaserDrillingAnalysis(model, parameters)
        simulation.Run()

        print(f"------------ Completed simulation for: {project_parameters_file_name} -------------")
