import h5py
import numpy as np

# Create an HDF5 file
filename = "./data/example.h5"
with h5py.File(filename, "w") as file:
    # Create a group
    group = file.create_group("my_group")

    # Create a dataset with some data
    data = np.array([[1, 2, 3], [4, 5, 6]])
    dataset = group.create_dataset("my_dataset", data=data)

    # Create an attribute
    attribute_value = "This is an attribute."
    group.attrs["my_attribute"] = attribute_value

    # You can also create additional groups, datasets, and attributes as needed.
    # ...

# Reading data from the created HDF5 file (optional)
with h5py.File(filename, "r") as file:
    # Access the dataset
    loaded_data = file["my_group/my_dataset"][:]
    print("Loaded Data:")
    print(loaded_data)

    # Access the attribute
    loaded_attribute = file["my_group"].attrs["my_attribute"]
    print("Loaded Attribute:")
    print(loaded_attribute)
