import scipy.io
import numpy as np

def mat_to_txt(mat_file, txt_file):
    # Load the .mat file
    mat_data = scipy.io.loadmat(mat_file)

    # Find the variable name you want to export
    var_name = None
    for key in mat_data.keys():
        if not key.startswith("__"):
            var_name = key
            break

    if var_name is None:
        print("No variable found in .mat file")
        return

    # Get the data
    data = mat_data[var_name]

    # Save the data to a .txt file
    np.savetxt(txt_file, data, fmt='%s', delimiter='\t')

# Replace these with the paths to your input .mat file and output .txt file
input_mat_file = "r1_tf1_new.mat"
output_txt_file = "r1_tf1_new.txt"

mat_to_txt(input_mat_file, output_txt_file)