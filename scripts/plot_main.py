import pandas as pd
import matplotlib.pyplot as plt

# Read data from the file
<<<<<<< HEAD
data = pd.read_csv("results/data.txt")
=======
data = pd.read_csv("/home/raghav/cws/temp/axbycz_probabilistic_method/results/plots_data.txt")
>>>>>>> main

# Create the plot
plt.figure()

<<<<<<< HEAD
plt.plot(data["Scramble Rate"], data["Error 1"], "o-r", label="Method 1")
plt.plot(data["Scramble Rate"], data["Error 3"], "s-b", label="Method 3")
=======
plt.plot(data["Scramble Rate"].to_numpy(), data["Error 1"].to_numpy(), "o-r", label="Prob 1")
plt.plot(data["Scramble Rate"].to_numpy(), data["Error 3"].to_numpy(), "s-b", label="Iterative")
>>>>>>> main

# Set labels and title
plt.xlabel("Scramble Rate (%)")
plt.ylabel("Error")
plt.title("Real Data")

<<<<<<< HEAD
# Add legend and grid
plt.legend()
plt.grid(True)

# Save the plot to a file
plt.savefig("results/Error_vs_Scramble_Rate_Python.png")
=======
# Add legend inside the plot box at the top right
plt.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95))

# Add grid
#plt.grid(True)

# Save the plot to a file
plt.savefig("/home/raghav/cws/temp/axbycz_probabilistic_method/results/Error_vs_Scramble_Rate_Python.png")
>>>>>>> main

# Show the plot
plt.show()