import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv("/home/raghav/cws/axbycz_probabilistic_method/results/results.txt", sep='\t')

plt.plot(data["Num_Matrices"].to_numpy(), data["Error1"].to_numpy(), "o-r", label="Probabilistic")
plt.plot(data["Num_Matrices"].to_numpy(), data["Error3"].to_numpy(), "s-b", label="Iterative")
plt.xlabel("Number of Data Triplets")
plt.ylabel("Error")
#plt.title("Error vs Number of Data Triplets")
plt.legend()
#plt.grid(True)
plt.savefig("/home/raghav/cws/axbycz_probabilistic_method/results/Error_vs_Datasets_3.png")
plt.show()