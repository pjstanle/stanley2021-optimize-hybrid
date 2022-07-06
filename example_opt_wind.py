import numpy as np
import gfopt.greedy_algorithm
import gfopt.genetic_algorithm

import openmdao.api as om

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
import scipy

import os
from dotenv import load_dotenv

from hybrid.sites import SiteInfo, flatirons_site
from hybrid.wind_source import WindPlant
from hybrid.keys import set_developer_nrel_gov_key

def get_turbine_locs(turbine_exists):

    global xlocs
    global ylocs

    layout_x = np.zeros(int(np.sum(turbine_exists)))
    layout_y = np.zeros(int(np.sum(turbine_exists)))
    ind = 0
    for i in range(len(turbine_exists)):
        if int(turbine_exists[i]) == 1:
            layout_x[ind] = xlocs[i]
            layout_y[ind] = ylocs[i]
            ind += 1

    return layout_x, layout_y


def calc_spacing(layout_x,layout_y):

        nTurbs = len(layout_x)
        npairs = int((nTurbs*(nTurbs-1))/2)
        spacing = np.zeros(npairs)

        ind = 0
        for i in range(nTurbs):
            for j in range(i,nTurbs):
                if i != j:
                    spacing[ind] = np.sqrt((layout_x[i]-layout_x[j])**2+(layout_y[i]-layout_y[j])**2)
                    ind += 1

        return spacing


def AEP_obj(x):

    global wind_plant
    global function_calls
    global side
    global minSpacing
    global rotor_diameter

    turbine_exists = x[:]

    layout_x, layout_y = get_turbine_locs(turbine_exists)

    if len(layout_x) > 1:
        spacing = calc_spacing(layout_x,layout_y)
        spacing_con = np.min(spacing) - minSpacing*rotor_diameter
    elif len(layout_x) == 1:
        spacing_con = 0.0
    else:
        spacing_con = -1.0

    if len(layout_x) < 1 or np.max(layout_x) > side or np.max(layout_y) > side or np.min(layout_x) < 0.0 or np.min(layout_y) < 0.0 or spacing_con < 0.0:
        return 1E6

    else:
        function_calls += 1
        wind_plant.modify_coordinates(layout_x,layout_y)
        wind_plant.simulate(1)
        # AEP = wind_plant.annual_energy_kw()/1E8    
        AEP = wind_plant.system_model.Outputs.annual_energy      

        if AEP/len(layout_x) < 4900000.0:
            return 1E6
        else:
            return -AEP


## Main
if __name__ == "__main__":

    global wind_plant
    global function_calls
    global turbine_exists
    global side
    global minSpacing
    global rotor_diameter

    global xlocs
    global ylocs

    # Set API key
    load_dotenv()
    NREL_API_KEY = os.getenv("NREL_API_KEY")
    set_developer_nrel_gov_key(NREL_API_KEY)  # Set this key manually here if you are not setting it using the .env

    # Get resource
    lat = flatirons_site['lat']
    lon = flatirons_site['lon']
    site = SiteInfo(flatirons_site)

    # Create model
    wind_plant = WindPlant(site, 0.0)
    wind_plant.wake_model = 1


    function_calls = 0

    side = 1200.0
    edges = np.array([0.0,0.0, 0.0,side, side,side, side,0.0]).reshape(4,2)
    grid_size = 10

    x_grid = np.linspace(0.0,side,grid_size)
    y_grid = np.linspace(0.0,side,grid_size)

    xlocs = np.zeros(grid_size*grid_size)
    ylocs = np.zeros(grid_size*grid_size)
    for i in range(grid_size):
        for j in range(grid_size):
            if i%2 == 0:
                xlocs[i*grid_size+j] = x_grid[j]
                ylocs[i*grid_size+j] = y_grid[i]
            else:
                xlocs[(i+1)*grid_size-j-1] = x_grid[j]
                ylocs[i*grid_size+j] = y_grid[i]

    minSpacing = 2.0 #rotor diameters
    rotor_diameter = 77.0

    # gradient-free optimization
    start_time = time.time()
    function_calls = 0

    ga = gfopt.greedy_algorithm.GreedyAlgorithm()
    # ga = gfopt.genetic_algorithm.GeneticAlgorithm()
    ga.bits = np.zeros(grid_size*grid_size,dtype=int)
    ga.bounds = np.zeros((grid_size*grid_size,2))
    ga.variable_type = np.array([])
    for i in range(grid_size*grid_size):
        ga.variable_type = np.append(ga.variable_type,"int")
        ga.bounds[i] = (0,1)
    ga.objective_function = AEP_obj

    # ga.population_size = 100
    # ga.convergence_iters = 25
    # ga.max_generation = 1000

    ga.optimize_switch(initialize=5)
    # ga.optimize_ga(initialize=5)

    opt_val = ga.optimized_function_value
    DVopt = ga.optimized_design_variables

    run_time = time.time() - start_time
    print("opt_val: ", opt_val)
    print("time to run: ", run_time)
    print("function calls: ", function_calls)

    

    xf, yf = get_turbine_locs(DVopt)    
    print("number of turbines: ", len(xf))

    plt.plot([0.0,side,side,0.0,0.0],[0.0,0.0,side,side,0.0],"--k",linewidth=0.5)
    ax = plt.gca()
    for i in range(len(xf)):
        circ = plt.Circle((xf[i],yf[i]),rotor_diameter/2.0,facecolor='C0',edgecolor='C0')
        ax.add_patch(circ)

    plt.axis("equal")
    plt.axis("off")
    plt.show()    


