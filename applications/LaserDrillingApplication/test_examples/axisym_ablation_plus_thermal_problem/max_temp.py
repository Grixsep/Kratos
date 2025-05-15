import h5py
import numpy as np
import matplotlib.pyplot as plt

# Open the HDF5 file (replace 'your_file.h5' with your actual file path)
file_path = "temperature_db.hdf5"

# Lists to store results
times = []
max_temperatures = []
max_temp_points = []

# Read the HDF5 file
with h5py.File(file_path, "r") as f:
    # Iterate through all top-level groups (time steps)
    for group_name in f.keys():
        group = f[group_name]
        if isinstance(group, h5py.Group):
            # Get the time from the group's attribute
            time_str = group.attrs.get("time")
            if time_str is None:
                print(f"Warning: No 'time' attribute in group {group_name}")
                continue
            try:
                time = float(time_str)
            except ValueError:
                print(f"Warning: Cannot convert 'time' attribute to float for group {group_name}")
                continue

            # Access the dataset inside the group (same name as group)
            dataset_name = group_name
            if dataset_name not in group:
                print(f"Warning: No dataset '{dataset_name}' in group {group_name}")
                continue
            dataset = group[dataset_name]

            # Check if it's a dataset
            if not isinstance(dataset, h5py.Dataset):
                print(f"Warning: '{dataset_name}' in group {group_name} is not a dataset")
                continue

            # Load the data (structured array)
            data = dataset[:]

            # Check if the dataset has the expected fields
            expected_fields = ["Id", "x_coord", "y_coord", "temperature"]
            if not all(field in data.dtype.names for field in expected_fields):
                print(f"Warning: Dataset in group {group_name} does not have the expected fields")
                continue

            # Extract temperature and Id fields
            temperatures = data["temperature"]
            ids = data["Id"]

            # Find the maximum temperature and corresponding point ID
            if len(temperatures) > 0:
                max_temp = np.max(temperatures)
                max_idx = np.argmax(temperatures)
                point_id = ids[max_idx]
            else:
                print(f"Warning: No data in dataset for group {group_name}")
                continue

            # Store results
            times.append(time)
            max_temperatures.append(max_temp)
            max_temp_points.append(point_id)
            print(f"Time {time}: Max temp = {max_temp:.2f} at point ID {int(point_id)}")

# Check if any data was found
if not times:
    print("Error: No valid datasets found in the HDF5 file.")
    exit(1)

# Sort by time
sorted_indices = np.argsort(times)
times = np.array(times)[sorted_indices]
max_temperatures = np.array(max_temperatures)[sorted_indices]

# Create the plot
plt.figure(figsize=(10, 6))
plt.plot(times, max_temperatures, ".", "b-", label="Max Temperature")
plt.xlabel("Time")
plt.ylabel("Temperature")
plt.title("Evolution of Mmximum temperature over time")
plt.grid(False)
plt.legend()

# Save the plot
# plt.savefig("max_temperature_plot.png")
# print("Plot saved as 'max_temperature_plot.png'")
plt.show()
plt.close()
