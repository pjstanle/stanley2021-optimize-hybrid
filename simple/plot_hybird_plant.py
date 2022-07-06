import numpy as np
from shapely.geometry import Polygon, Point
import matplotlib.pyplot as plt

def plot_hybrid(x):
    global solar_width
    global solar_height

    plt.cla()
    turbine_x, turbine_y, solar_x, solar_y = create_farm(x)
    for i in range(len(turbine_x)):
        turb = plt.Circle((turbine_x[i],turbine_y[i]),59.0,facecolor="C0")
        plt.gca().add_patch(turb)

    for i in range(len(solar_x)):
        panel = plt.Rectangle((solar_x[i]-solar_width/2.0,solar_y[i]-solar_height/2.0),solar_width,solar_height,facecolor="C1")
        plt.gca().add_patch(panel)
    plt.axis("equal")
    plt.pause(0.0001)


def create_farm(plant_array):

    global xlocs
    global ylocs

    nturbs = np.count_nonzero(plant_array==1)
    narrays = np.count_nonzero(plant_array==2)
    turbine_x = np.zeros(nturbs)
    turbine_y = np.zeros(nturbs)
    solar_x = np.zeros(narrays)
    solar_y = np.zeros(narrays)

    turbine_ind = 0
    solar_ind = 0
    for i in range(len(plant_array)):
        if int(plant_array[i]) == 1:
            turbine_x[turbine_ind] = xlocs[i]
            turbine_y[turbine_ind] = ylocs[i]
            turbine_ind += 1

        if int(plant_array[i]) == 2:
            solar_x[solar_ind] = xlocs[i]
            solar_y[solar_ind] = ylocs[i]
            solar_ind += 1

    return turbine_x, turbine_y, solar_x, solar_y


if __name__=="__main__":

    global solar_width
    global solar_height
    global xlocs
    global ylocs

    side = 1600.0
    grid_size = 20
    x_grid = np.linspace(0.0,side,grid_size)
    y_grid = np.linspace(0.0,side,grid_size)

    xlocs, ylocs = np.meshgrid(x_grid,y_grid)
    xlocs = np.ndarray.flatten(xlocs)
    ylocs = np.ndarray.flatten(ylocs)

    solar_width = x_grid[1]-x_grid[0]-1.0
    solar_height = y_grid[1]-y_grid[0]-1.0

    x = np.array([0.000, 0.000, 2.000, 2.000, 0.000, 2.000, 0.000, 2.000, 2.000,
       2.000, 0.000, 0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000,
       0.000, 0.000, 2.000, 0.000, 0.000, 2.000, 2.000, 2.000, 0.000,
       2.000, 0.000, 2.000, 2.000, 0.000, 0.000, 0.000, 0.000, 0.000,
       0.000, 2.000, 0.000, 2.000, 0.000, 2.000, 0.000, 0.000, 0.000,
       2.000, 0.000, 0.000, 2.000, 2.000, 2.000, 0.000, 2.000, 0.000,
       0.000, 0.000, 0.000, 2.000, 2.000, 2.000, 2.000, 0.000, 0.000,
       0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.000, 2.000, 2.000,
       0.000, 2.000, 0.000, 2.000, 0.000, 2.000, 0.000, 0.000, 0.000,
       0.000, 0.000, 1.000, 0.000, 0.000, 0.000, 0.000, 2.000, 0.000,
       2.000, 2.000, 2.000, 2.000, 2.000, 0.000, 0.000, 0.000, 0.000,
       0.000, 2.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.000, 2.000,
       2.000, 2.000, 2.000, 2.000, 2.000, 0.000, 0.000, 0.000, 0.000,
       0.000, 0.000, 1.000, 2.000, 2.000, 0.000, 0.000, 0.000, 2.000,
       0.000, 0.000, 0.000, 0.000, 2.000, 2.000, 0.000, 2.000, 2.000,
       2.000, 2.000, 0.000, 0.000, 0.000, 2.000, 2.000, 0.000, 2.000,
       2.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.000, 0.000,
       0.000, 2.000, 2.000, 2.000, 2.000, 0.000, 0.000, 2.000, 2.000,
       0.000, 2.000, 0.000, 0.000, 0.000, 0.000, 1.000, 0.000, 0.000,
       0.000, 0.000, 2.000, 2.000, 2.000, 2.000, 2.000, 0.000, 0.000,
       2.000, 2.000, 2.000, 2.000, 2.000, 2.000, 0.000, 0.000, 0.000,
       0.000, 0.000, 2.000, 0.000, 2.000, 2.000, 0.000, 0.000, 2.000,
       2.000, 0.000, 2.000, 0.000, 0.000, 0.000, 2.000, 2.000, 2.000,
       0.000, 0.000, 0.000, 0.000, 0.000, 2.000, 0.000, 2.000, 2.000,
       2.000, 2.000, 2.000, 2.000, 0.000, 0.000, 0.000, 0.000, 0.000,
       2.000, 2.000, 0.000, 2.000, 2.000, 2.000, 0.000, 2.000, 2.000,
       2.000, 2.000, 2.000, 2.000, 2.000, 0.000, 0.000, 0.000, 1.000,
       0.000, 0.000, 0.000, 0.000, 0.000, 2.000, 2.000, 2.000, 2.000,
       0.000, 0.000, 0.000, 2.000, 2.000, 0.000, 2.000, 2.000, 0.000,
       0.000, 0.000, 0.000, 0.000, 0.000, 2.000, 2.000, 2.000, 0.000,
       2.000, 0.000, 2.000, 2.000, 0.000, 2.000, 2.000, 0.000, 2.000,
       2.000, 2.000, 0.000, 0.000, 0.000, 2.000, 2.000, 2.000, 2.000,
       2.000, 0.000, 0.000, 2.000, 0.000, 2.000, 2.000, 2.000, 2.000,
       2.000, 2.000, 2.000, 2.000, 2.000, 2.000, 0.000, 2.000, 2.000,
       2.000, 2.000, 2.000, 2.000, 2.000, 0.000, 2.000, 2.000, 2.000,
       0.000, 0.000, 0.000, 0.000, 2.000, 2.000, 2.000, 0.000, 0.000,
       0.000, 0.000, 2.000, 0.000, 2.000, 2.000, 0.000, 2.000, 0.000,
       0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.000, 0.000, 0.000,
       0.000, 0.000, 0.000, 0.000, 2.000, 0.000, 2.000, 0.000, 0.000,
       0.000, 0.000, 0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000,
       0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000,
       0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
       0.000, 0.000, 2.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
       0.000, 0.000, 0.000, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000,
       0.000, 0.000, 2.000, 0.000])

    plot_hybrid(x)
    plt.show()