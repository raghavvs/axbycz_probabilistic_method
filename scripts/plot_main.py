import pandas as pd
import matplotlib.pyplot as plt

# Read data from the file
data = pd.read_csv("results/data.txt")

# Create the plot
plt.figure()

plt.plot(data["Scramble Rate"], data["Error 1"], "o-r", label="Method 1")
plt.plot(data["Scramble Rate"], data["Error 3"], "s-b", label="Method 3")

# Set labels and title
plt.xlabel("Scramble Rate (%)")
plt.ylabel("Error")
plt.title("Real Data")

# Add legend and grid
plt.legend()
plt.grid(True)

# Save the plot to a file
plt.savefig("results/Error_vs_Scramble_Rate_Python.png")

# Show the plot
plt.show()