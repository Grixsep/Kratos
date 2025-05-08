import json
import os

# List of beam_waist_diameter values
beam_waist_diameters = [0.02, 0.07, 0.12]

# Load the input JSON file
with open("ProjectParameters.json", "r") as file:
    data = json.load(file)

# Create output directory if it doesn't exist
output_dir = "parameter_jsons"
os.makedirs(output_dir, exist_ok=True)

# Generate a new JSON file for each beam_waist_diameter value
for i, diameter in enumerate(beam_waist_diameters):
    # Create a deep copy of the original data
    modified_data = json.loads(json.dumps(data))

    # Modify the beam_waist_diameter value
    modified_data["problem_data"]["beam_waist_diameter"] = diameter

    # Modify the output_name in gid_output
    modified_data["output_processes"]["gid_output"][0]["Parameters"]["output_name"] = (
        f"gid_output/waist_{i}/axisym_ablation_plus_thermal_problem"
    )

    # Write the modified data to a new JSON file
    output_filename = os.path.join(output_dir, f"ProjectParameters_{i}.json")
    with open(output_filename, "w") as file:
        json.dump(modified_data, file, indent=4)

    print(
        f"Created {output_filename} with beam_waist_diameter = {diameter} and output_name = gid_output/waist_{i}/axisym_ablation_plus_thermal_problem"
    )
