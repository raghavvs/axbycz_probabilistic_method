import pandas as pd
import matplotlib.pyplot as plt

# Read data from the file
data = pd.read_csv("/home/raghav/cws/temp/axbycz_probabilistic_method/results/plots_data.txt")

# Create the plot
plt.figure()

plt.plot(data["Scramble Rate"].to_numpy(), data["Error 1"].to_numpy(), "o-r", label="Prob 1")
plt.plot(data["Scramble Rate"].to_numpy(), data["Error 3"].to_numpy(), "s-b", label="Iterative")

# Set labels and title
plt.xlabel("Scramble Rate (%)")
plt.ylabel("Error")
plt.title("Real Data")

# Add legend inside the plot box at the top right
plt.legend(loc='upper right', bbox_to_anchor=(0.95, 0.95))

# Add grid
#plt.grid(True)

# Save the plot to a file
plt.savefig("/home/raghav/cws/temp/axbycz_probabilistic_method/results/Error_vs_Scramble_Rate_Python.png")

# Show the plot
plt.show()