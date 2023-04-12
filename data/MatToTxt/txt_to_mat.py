import numpy as np
import scipy.io

def txt_to_mat(txt_file, mat_file):
    # Load the .txt file
    data = np.loadtxt(txt_file, delimiter=' ')

    # Save the data to a .mat file
    scipy.io.savemat(mat_file, {'data': data})

# Replace these with the paths to your input .txt file and output .mat file
input_txt_file = "A_abc.txt"
output_mat_file = "A_abc.mat"

txt_to_mat(input_txt_file, output_mat_file)