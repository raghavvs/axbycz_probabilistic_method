import matplotlib.pyplot as plt

def plot_boxplot(data, points):
    fig, ax = plt.subplots()
    ax.boxplot(data, positions=points)
    plt.show()