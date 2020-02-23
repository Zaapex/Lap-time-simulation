import pandas as pd
import webbrowser
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def show_data(name_of_track):
        """Present data in beautiful manner in my browser"""

        df = pd.read_csv(name_of_track, index_col=0, low_memory=False)

        f = open(name_of_track + ".html", "x")
        f.write(df.to_html())  # create an HTML file of your data set
        webbrowser.open(name_of_track+".html")

def show_track(name_of_track):
        """Parametric plot of x and y values-> you can see the track"""

        df = pd.read_csv(name_of_track, index_col=0)

        x = df["x"][:]  # simple x and y position of the track
        y = df["y"][:]

        plt.plot(x[:2500], y[:2500], "r.")
        plt.plot(x[2501:5000], y[2501:5000], "g.")
        plt.plot(x[5001:7500], y[5001:7500], "b.")
        plt.plot(x[7501:10000], y[7501:10000], "y.")
        plt.show()

def tire_data(tire_path):
        """Collect data from tire and return value that are interesting, Podatki\Pnevmatike\Tire_Hoosier.csv"""

        df = pd.read_csv(tire_path, index_col=0, sep=";", low_memory=False)

        return df

def acc_or_brake(name_of_track):
        """Show on the map where we accelerate and where we are braking"""

        df = name_of_track

        x = df["x"]  # x, y track data and corresponding acceleration in that point
        y = df["y"]
        a = df["acceleration"]

        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        fig, axs = plt.subplots(2, 1, sharex=True, sharey=True)

        mean = np.mean(a)  # defining the median and standard deviation of the acceleration
        std = np.std(a)

        norm = plt.Normalize(mean - std, mean + std)
        lc = LineCollection(segments, cmap='viridis', norm=norm)

        lc.set_array(a)
        lc.set_linewidth(2)
        line = axs[0].add_collection(lc)
        fig.colorbar(line, ax=axs[0])

        cmap = ListedColormap(['r', 'g', 'b'])
        norm = BoundaryNorm([-10, -1, 1, 10], cmap.N)  # your set of boundaries
        lc = LineCollection(segments, cmap=cmap, norm=norm)
        lc.set_array(a)
        lc.set_linewidth(2)
        line = axs[1].add_collection(lc)
        fig.colorbar(line, ax=axs[1])

        axs[0].set_xlim(-42, 60)
        axs[0].set_ylim(-150, 5)
        plt.xlabel("x [m]")
        plt.ylabel("y [m]")
        plt.show()

def plot_3D(name_of_track, values):

    df = name_of_track

    a = values[0]  # usally x and y data
    b = values[1]
    c = values[2]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot(df[a], df[b], df[c], "r.")
    plt.show()


name = "FSG_Layout"
tires_data = "Podatki\Pnevmatike\Tire_data.csv"
tire_name = "Tire_Hoosier"



