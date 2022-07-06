import numpy as np
from shapely.geometry import Polygon, Point, MultiPolygon
from gradient_free import GeneticAlgorithm, GreedyAlgorithm
from model_hybrid import SimpleHybrid
import floris.tools as wfct
import matplotlib.pyplot as plt
from hybrid_optimization_tools import hybrid_sweep, hybrid_sweep_w_battery
import time
import geopandas as gpd
from shapely.ops import unary_union

def create_farm_stamp(plant_array):

    global xlocs
    global ylocs

    global stamp_x
    global stamp_y

    nstamps = len(stamp_x)
    nturbs = np.count_nonzero(plant_array==1)*nstamps
    narrays = np.count_nonzero(plant_array==2)*nstamps
    turbine_x = np.zeros(nturbs)
    turbine_y = np.zeros(nturbs)
    solar_x = np.zeros(narrays)
    solar_y = np.zeros(narrays)

    turbine_ind = 0
    solar_ind = 0
    for k in range(nstamps):
        for i in range(len(plant_array)):
            if int(plant_array[i]) == 1:
                turbine_x[turbine_ind] = xlocs[i] + stamp_x[k]
                turbine_y[turbine_ind] = ylocs[i] + stamp_y[k]
                turbine_ind += 1

            if int(plant_array[i]) == 2:
                solar_x[solar_ind] = xlocs[i] + stamp_x[k]
                solar_y[solar_ind] = ylocs[i] + stamp_y[k]
                solar_ind += 1

    return turbine_x, turbine_y, solar_x, solar_y


def plot_hybrid(turbine_x,turbine_y,solar_x,solar_y):
    global solar_width
    global solar_height

    plt.cla()
    for i in range(len(turbine_x)):
        turb = plt.Circle((turbine_x[i],turbine_y[i]),59.0,facecolor="C0")
        plt.gca().add_patch(turb)

    for i in range(len(solar_x)):
        panel = plt.Rectangle((solar_x[i]-solar_width/2.0,solar_y[i]-solar_height/2.0),solar_width,solar_height,facecolor="C1")
        plt.gca().add_patch(panel)
    plt.axis("equal")
    plt.pause(0.0001)


def union_solar(solar_x,solar_y):
    global solar_width
    global solar_height

    narrays = len(solar_x)
    panels = [0]*narrays
    for i in range(narrays):
        panels[i] = Polygon(([solar_x[i]-solar_width/2.0,solar_y[i]-solar_height/2.0],[solar_x[i]+solar_width/2.0,solar_y[i]-solar_height/2.0],[solar_x[i]+solar_width/2.0,solar_y[i]+solar_height/2.0],[solar_x[i]-solar_width/2.0,solar_y[i]+solar_height/2.0]))

    panel_union = unary_union(panels)

    ngroups = len(panel_union)
    areas = np.zeros(ngroups)
    for i in range(ngroups):
        areas[i] = panel_union[i].area
    
    return ngroups, areas


if __name__=="__main__":

    global xlocs
    global ylocs
    global solar_width
    global solar_height
    global stamp_x
    global stamp_y

    grid_size = 5
    side = 118.0*(grid_size-1)
    x_grid = np.linspace(0.0,side,grid_size)
    y_grid = np.linspace(0.0,side,grid_size)

    solar_width = x_grid[1]-x_grid[0]
    solar_height = y_grid[1]-y_grid[0]

    xlocs, ylocs = np.meshgrid(x_grid,y_grid)
    xlocs = np.ndarray.flatten(xlocs)
    ylocs = np.ndarray.flatten(ylocs)

    farm_width = max(xlocs)-min(xlocs)
    farm_height = max(ylocs)-min(ylocs)
    dx = x_grid[1]-x_grid[0]
    dy = y_grid[1]-y_grid[0]

    stx = np.array([0.,1.,2.])*(farm_width+dx)
    sty = np.array([0.,1.,2.])*(farm_height+dy)
    stamp_x,stamp_y = np.meshgrid(stx,sty)
    stamp_x = np.ndarray.flatten(stamp_x)
    stamp_y = np.ndarray.flatten(stamp_y)


    nlocs = len(xlocs)
    np.random.seed(3)
    plant_array = np.random.randint(0,high=3,size=nlocs)
    c = 0
    for i in range(nlocs):
        if plant_array[i]==1:
            if c != 2:
                plant_array[i]=0
            c = (c+1)%3


    
    turbine_x, turbine_y, solar_x, solar_y = create_farm_stamp(plant_array)

    plot_hybrid(turbine_x, turbine_y, solar_x, solar_y)

    ngroups, areas = union_solar(solar_x,solar_y)
    print(ngroups,-np.sort(-areas))


    plt.show()
