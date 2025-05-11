import os
import sys
import json
import KratosMultiphysics
from KratosMultiphysics.LaserDrillingApplication.laserdrilling_analysis import LaserDrillingAnalysis


def run_simulation(project_parameters_file, subdirectory):
    """Run a laser drilling simulation using the given ProjectParameters.json file."""
    try:
        # Read the ProjectParameters.json file
        with open(project_parameters_file, "r") as file:
            parameters = KratosMultiphysics.Parameters(file.read())

        # Create a new model
        model = KratosMultiphysics.Model()

        # Create and run the LaserDrillingAnalysis
        simulation = LaserDrillingAnalysis(model, parameters)

        print(f"[LaserDrillingSimulation] Starting simulation in {subdirectory}")
        simulation.Run()
        print(f"[LaserDrillingSimulation] Simulation in {subdirectory} completed successfully")

    except FileNotFoundError:
        print(f"[LaserDrillingSimulation] ERROR: ProjectParameters.json not found in {subdirectory}")
    except json.JSONDecodeError:
        print(f"[LaserDrillingSimulation] ERROR: Invalid JSON in ProjectParameters.json in {subdirectory}")
    except KeyError as e:
        print(
            f"[LaserDrillingSimulation] ERROR: Missing parameter {str(e)} in ProjectParameters.json in {subdirectory}"
        )
    except Exception as e:
        print(f"[LaserDrillingSimulation] ERROR: Error running simulation in {subdirectory}: {str(e)}")


def main(directory_path):
    """Main function to iterate through subdirectories and launch simulations."""
    # Check if the directory exists
    if not os.path.isdir(directory_path):
        print(f"[LaserDrillingSimulation] ERROR: Directory {directory_path} does not exist")
        sys.exit(1)

    # Iterate through subdirectories
    for subdir in os.listdir(directory_path):
        subdir_path = os.path.join(directory_path, subdir)
        if os.path.isdir(subdir_path):
            project_parameters_file = os.path.join(subdir_path, "ProjectParameters.json")
            if os.path.isfile(project_parameters_file):
                print(f"[LaserDrillingSimulation] Processing subdirectory: {subdir}")
                run_simulation(project_parameters_file, subdir)
            else:
                print(f"[LaserDrillingSimulation] WARNING: No ProjectParameters.json found in {subdir}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("[LaserDrillingSimulation] ERROR: Usage: python script.py <directory_path>")
        sys.exit(1)

    directory_path = sys.argv[1]
    main(directory_path)
