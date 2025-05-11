import json
import os
import sys

# List of beam_waist_diameter values
beam_waist_diameters = [0.02, 0.071, 0.12]

# Load the input JSON files
project_parameters_input_filename = os.path.join("parameters", "ProjectParameters.json")
laser_parameters_input_filename = os.path.join("parameters", "LaserParameters.json")
materials_parameters_input_filename = os.path.join("parameters", "MaterialsParameters.json")
environment_parameters_input_filename = os.path.join("parameters", "EnvironmentParameters.json")

with open(project_parameters_input_filename, "r") as input_project_file:
    project_data = json.load(input_project_file)
# Create a deep copy of the original data
modified_project_data = json.loads(json.dumps(project_data))

with open(laser_parameters_input_filename, "r") as input_laser_file:
    laser_data = json.load(input_laser_file)
# Create a deep copy of the original data
modified_laser_data = json.loads(json.dumps(laser_data))

with open(materials_parameters_input_filename, "r") as input_materials_file:
    materials_data = json.load(input_materials_file)
# Create a deep copy of the original data
modified_materials_data = json.loads(json.dumps(materials_data))

with open(environment_parameters_input_filename, "r") as input_environment_file:
    environment_data = json.load(input_environment_file)
# Create a deep copy of the original data
modified_environment_data = json.loads(json.dumps(environment_data))


# Create output directory if it doesn't exist
batches_output_dir = "waist_varying_batch"
try:
    os.makedirs(batches_output_dir, exist_ok=False)
except FileExistsError:
    print("The directory for the batch run already exists")
    sys.exit(1)


# Generate new parameters files for each beam_waist_diameter value
for i, beam_waist_diam in enumerate(beam_waist_diameters):
    # Create an output directory for each parameter configuration
    batch_dir = os.path.join(batches_output_dir, f"beam_waist_{i}")
    os.makedirs(batch_dir, exist_ok=True)

    # Modify ProjectParameters
    # Modify the output_name in gid_output
    modified_project_data["output_processes"]["gid_output"][0]["Parameters"]["output_name"] = os.path.join(
        batch_dir, "gid_output", "axisym_ablation_plus_thermal_problem"
    )

    # Modify the path of the parameter files
    modified_project_data["solver_settings"]["material_import_settings"]["materials_filename"] = os.path.join(
        batch_dir, "MaterialsParameters.json"
    )
    modified_project_data["solver_settings"]["material_import_settings"]["laser_filename"] = os.path.join(
        batch_dir, "LaserParameters.json"
    )
    modified_project_data["solver_settings"]["material_import_settings"]["environment_filename"] = os.path.join(
        batch_dir, "EnvironmentParameters.json"
    )

    # Output the new ProjectParameters
    output_project_filename = os.path.join(batch_dir, "ProjectParameters.json")
    with open(output_project_filename, "w") as output_project_file:
        json.dump(modified_project_data, output_project_file, indent=4)

    # Modify LaserParameters
    # Modify the beam_waist_diameter value
    modified_laser_data["properties"][0]["Material"]["Variables"]["beam_waist_diameter"] = beam_waist_diam

    # Output the new LaserParameters
    output_laser_filename = os.path.join(batch_dir, "LaserParameters.json")
    with open(output_laser_filename, "w") as output_laser_file:
        json.dump(modified_laser_data, output_laser_file, indent=4)

    # Modify MaterialsParameters
    # Do nothing in this case
    # Output the new MaterialsParameters
    output_materials_filename = os.path.join(batch_dir, "MaterialsParameters.json")
    with open(output_materials_filename, "w") as output_materials_file:
        json.dump(modified_materials_data, output_materials_file, indent=4)

    # Modify EnvironmentParameters
    # Do nothing in this case
    # Output the new EnvironmentParameters
    output_environment_filename = os.path.join(batch_dir, "EnvironmentParameters.json")
    with open(output_environment_filename, "w") as output_environment_file:
        json.dump(modified_environment_data, output_environment_file, indent=4)
