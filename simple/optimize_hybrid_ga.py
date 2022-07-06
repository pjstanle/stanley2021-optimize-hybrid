import numpy as np
from shapely.geometry import Polygon, Point
from gradient_free import GeneticAlgorithm, GreedyAlgorithm
from model_hybrid import SimpleHybrid
import floris.tools as wfct
import matplotlib.pyplot as plt
from hybrid_optimization_tools import hybrid_greedy
import time


def initialize_plant(npoints,wind_probability,solar_probability):

    plant_array = np.zeros(npoints,dtype=int)
    rand_array = np.random.rand(npoints)
    for i in range(npoints):
        if rand_array[i] < wind_probability:
            plant_array[i] = 1
        elif rand_array[i] < wind_probability+solar_probability:
            plant_array[i] = 2
    return plant_array


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


def AEP_obj(x):

    global plant
    global wind_wind_spacing
    global solar_solar_spacing
    global wind_solar_spacing

    global xlocs
    global ylocs

    global solar_width
    global solar_height

    plant_array = x[:]
    # solar_width = x[-2]
    # solar_height = x[-1]

    turbine_x, turbine_y, solar_x, solar_y = create_farm(plant_array)

    plant.turbine_x = turbine_x
    plant.turbine_y = turbine_y
    plant.solar_x = solar_x
    plant.solar_y = solar_y

    width = np.ones(len(solar_x))*solar_width
    height = np.ones(len(solar_x))*solar_height
    plant.solar_width = width
    plant.solar_height = height

    wind_solar = plant.get_wind_solar_distance()
    wind_wind = plant.get_wind_wind_distance()
    solar_solar = plant.get_solar_solar_distance()
    if len(wind_solar) == 0:
        wind_solar = np.array([wind_solar_spacing])
    if len(wind_wind) == 0:
        wind_wind = np.array([wind_wind_spacing])
    if len(solar_solar) == 0:
        solar_solar = np.array([solar_solar_spacing])

    if len(solar_x) == 0:
        solar_x = np.array([np.mean(xlocs)])
    if len(solar_y) == 0:
        solar_y = np.array([np.mean(ylocs)])

    boundary = Polygon(([min(xlocs)-solar_width/2.0,min(ylocs)-solar_height/2.0],[max(xlocs)+solar_width/2.0,min(ylocs)-solar_height/2.0],[max(xlocs)+solar_width/2.0,max(ylocs)+solar_height/2.0],[min(xlocs)-solar_width/2.0,max(ylocs)+solar_height/2.0]))
    bc = True
    for i in range(len(solar_x)):
        sa = Polygon(([solar_x[i]-solar_width/2.0,solar_y[i]-solar_height/2.0],[solar_x[i]+solar_width/2.0,solar_y[i]-solar_height/2.0],\
                    [solar_x[i]+solar_width/2.0,solar_y[i]+solar_height/2.0],[solar_x[i]-solar_width/2.0,solar_y[i]+solar_height/2.0]))
        if boundary.contains(sa)==False:
            bc = False

    if min(wind_solar) < wind_solar_spacing or min(wind_wind) < wind_wind_spacing or min(solar_solar) < solar_solar_spacing or bc == False: 
                # (max(solar_x)+solar_width/2.0) > max(xlocs) or (min(solar_x)-solar_width/2.0) < min(xlocs) or (max(solar_y)+solar_height/2.0) > max(ylocs) \
                # or (min(solar_y)-solar_height/2.0) < max(ylocs):
        return 1E20
    else:
        plant.get_plant_aep()
    return -plant.plant_aep


def capacity_obj(x):

    global plant
    global wind_wind_spacing
    global solar_solar_spacing
    global wind_solar_spacing

    global xlocs
    global ylocs

    global solar_width
    global solar_height

    global wd_array
    global ws_array
    global sr_array

    global interconnect

    global battery_storage
    global battery_output
    global min_energy

    plant_array = x[:]

    turbine_x, turbine_y, solar_x, solar_y = create_farm(plant_array)

    plant.turbine_x = turbine_x
    plant.turbine_y = turbine_y
    plant.solar_x = solar_x
    plant.solar_y = solar_y

    width = np.ones(len(solar_x))*solar_width
    height = np.ones(len(solar_x))*solar_height
    plant.solar_width = width
    plant.solar_height = height

    wind_solar = plant.get_wind_solar_distance()
    wind_wind = plant.get_wind_wind_distance()
    solar_solar = plant.get_solar_solar_distance()

    if (len(solar_x)==0 and len(turbine_x)==0):
        return 0.0
    elif len(solar_x)==0 or len(turbine_x)==0:
        if len(solar_x)==0:
            if min(wind_wind) < wind_wind_spacing:
                return 1E20
        elif len(turbine_x)==0:
            if min(solar_solar) < solar_solar_spacing:
                return 1E20
    else:
        if min(wind_solar) < wind_solar_spacing or min(wind_wind) < wind_wind_spacing or min(solar_solar) < solar_solar_spacing:
            return 1E20
        
    plant.get_solar_capacity()
    plant_capacity = sum(plant.solar_capacity_mw) + len(turbine_x)*2.5
    if plant_capacity > 10.0*interconnect:
        return 1E20
    else:
        nhours = len(wd_array)
        gross_energy = np.zeros(nhours)
        real_energy = np.zeros(nhours)
        for i in range(nhours):
            gross_energy[i] = plant.get_plant_power(wd_array[i],ws_array[i],sr_array[i])


            if gross_energy[i] > interconnect:
                real_energy[i] = interconnect
            else:
                real_energy[i] = gross_energy[i]
        
        curtailment = gross_energy-real_energy
        np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})

        battery_charge = np.zeros(nhours+1)
        real_energy_with_battery = np.zeros(nhours)
        for i in range(nhours):
            # charge the battery
            if curtailment[i] > 0.0 and battery_charge[i] < battery_storage:
                battery_charge[i+1] = battery_charge[i]+curtailment[i]
                real_energy_with_battery[i] = real_energy[i]
                if battery_charge[i+1] > battery_storage:
                    battery_charge[i+1] = battery_storage
            # discharge the battery
            elif real_energy[i] < min_energy:
                required_energy = min_energy-real_energy[i]
                provided_energy = min(required_energy,battery_charge[i])
                real_energy_with_battery[i] = real_energy[i]+provided_energy
                battery_charge[i+1] = battery_charge[i]-provided_energy
            else:
                battery_charge[i+1] = battery_charge[i]
                real_energy_with_battery[i] = real_energy[i]
            
        print("gross energy: ", gross_energy)
        print("real energy: ", real_energy)
        print("real energy sum: ", sum(real_energy))
        print("curtailment: ", curtailment)
        print("battery charge: ", battery_charge)
        print("real energy with battery: ", real_energy_with_battery)
        print("real energy with battery sum: ", sum(real_energy_with_battery))

        plt.figure(1)
        plt.plot(gross_energy,label="gross")
        plt.plot(real_energy,label="real")
        plt.plot(curtailment,label="curtail")
        plt.plot(battery_charge,label="charge")
        plt.plot(real_energy_with_battery,label="w/ battery")
        plt.legend()

        plt.figure(2)
        plt.plot(ws_array/12.0,label="wind")
        plt.plot(sr_array,label="solar")
        plt.legend()

        wo_bat = np.zeros(nhours)  
        w_bat = np.zeros(nhours)
        for i in range(nhours):
            wo_bat[i] = sum(real_energy[0:i])
            w_bat[i] = sum(real_energy_with_battery[0:i])
        
        plt.figure(3)
        plt.plot(wo_bat,label="w/o")
        plt.plot(w_bat,label="w/")

        plt.show()
        
        # if sum(curtailment) > 24.0:
        #     return 1E20
        # else:

        return -np.sum(real_energy_with_battery)
        


if __name__=="__main__":

    global xlocs
    global ylocs
    global plant
    global wind_wind_spacing
    global solar_solar_spacing
    global wind_solar_spacing

    global solar_width
    global solar_height

    global wd_array
    global ws_array
    global sr_array

    global interconnect

    global battery_storage
    global battery_output

    battery_storage = 3.0 # MWh
    battery_output = 3.0 # MW
    min_energy = 8.0

    interconnect = 10.0

    # wd_array = np.array([270.0])
    # ws_array = np.array([8.0])
    # sr_array = np.array([1.0])

    # nhours = 24
    # ws_array = np.random.rand(nhours)*10.0
    # sr_array = np.random.rand(nhours)
    ws_array = np.array([0.95,1.00,1.00,0.95,0.85,0.7,0.5,0.3,0.15,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.15,0.3,0.5,0.7,0.85])*6.0+6.0
    sr_array = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.15,0.3,0.5,0.7,0.85,0.95,1.00,1.00,0.95,0.85,0.7,0.5,0.3,0.15,0.05,0.0,0.0])
    ws_array = np.append(ws_array,ws_array)
    sr_array = np.append(sr_array,sr_array)
    # ws_array[35:37] = 0.0
    # sr_array[35:37] = 0.0
    wd_array = np.ones(len(ws_array))*270.0

    side = 800.0
    grid_size = 5
    x_grid = np.linspace(0.0,side,grid_size)
    y_grid = np.linspace(0.0,side,grid_size)

    # side = 800.0
    # grid_size = 5
    # x_grid = np.linspace(0.0,side,grid_size)
    # y_grid = np.linspace(0.0,side,grid_size)

    solar_width = x_grid[1]-x_grid[0]-1.0
    solar_height = y_grid[1]-y_grid[0]-1.0

    print(solar_width)

    xlocs, ylocs = np.meshgrid(x_grid,y_grid)
    xlocs = np.ndarray.flatten(xlocs)
    ylocs = np.ndarray.flatten(ylocs)
    
    plant = SimpleHybrid()
    plant.floris_model = wfct.floris_interface.FlorisInterface("/Users/astanley/Data/turbines/2_5mw_118d_88h.json")
    plant.floris_model.set_gch(False)

    plant.solar_capacity_factor = 0.2
    # plant.solar_capacity_factor = 0.20375836750705578
    plant.acres_per_MW = 5.0
    plant.shadow_x = 2.0 # rotor diameters of the sigma of shadow loss in x (east west)
    plant.shadow_y = 0.5 # rotor diameters of the sigma of shadow loss in y (north south)
    plant.shadow_scale = 0.5

    # plant.wind_directions = np.array([45.0])
    # plant.wind_speeds = np.array([8.0])
    # plant.wind_frequencies = np.array([0.5])
        
    D = 118.0
    wind_wind_spacing = D*2.0
    solar_solar_spacing = 0.0
    wind_solar_spacing = D*1.0



    nlocs = len(xlocs)
    ntech = 2

    # gradient-free optimization
    start_time = time.time()

    ga = GeneticAlgorithm()

    nvars = len(xlocs)
    ga.bits = np.zeros(nvars,dtype=int)
    ga.bounds = np.zeros((nvars,2))
    ga.variable_type = np.array([])
    for i in range(grid_size*grid_size):
        ga.bits[i] = 2
        ga.variable_type = np.append(ga.variable_type,"int")
        ga.bounds[i] = (0,3)
    ga.objective_function = capacity_obj

    ga.population_size = 50
    ga.convergence_iters = 10
    ga.max_generation = 1000

    # ga.optimize_ga(initialize="limit",crossover="chunk")
    # ga.optimize_ga(initialize="limit")

    # opt_val = ga.optimized_function_value
    # DVopt = ga.optimized_design_variables

    DVopt = np.array([2.000, 0.000, 5.000, 1.000, 4.000, 2.000, 0.000, 0.000, 0.000,
       1.000, 2.000, 2.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
       1.000, 0.000, 2.000, 3.000, 0.000, 4.000, 1.000])
    # DVopt = np.array([1., 0., 0., 2., 5., 0., 0., 0., 0., 1., 2., 0., 2., 2., 0., 0., 1.,
    #    0., 0., 2., 2., 0., 0., 1., 0.])
    # DVopt = np.array([2.000, 5.000, 1.000, 0.000, 1.000, 2.000, 2.000, 0.000, 1.000,
    #    0.000, 2.000, 2.000, 2.000, 0.000, 2.000, 2.000, 2.000, 4.000,
    #    2.000, 2.000, 2.000, 0.000, 1.000, 0.000, 2.000])
    opt_val = capacity_obj(DVopt)

    run_time = time.time() - start_time
    print("optimal design variables: ", repr(DVopt))
    print("opt_val: ", opt_val)
    print("time to run: ", run_time)

    turbine_x, turbine_y, solar_x, solar_y = create_farm(DVopt)
    for i in range(len(turbine_x)):
        turb = plt.Circle((turbine_x[i],turbine_y[i]),59.0,facecolor="C0")
        plt.gca().add_patch(turb)

    for i in range(len(solar_x)):
        dx = solar_width/2.0
        dy = solar_height/2.0
        cx = solar_x[i]
        cy = solar_y[i]
        x = [cx-dx,cx+dx,cx+dx,cx-dx,cx-dx]
        y = [cy-dy,cy-dy,cy+dy,cy+dy,cy-dy]
        plt.plot(x,y,color="C1",linewidth=2)
    
    plt.axis("equal")
    plt.show()


    # xf, yf = get_grid_locs(DVopt[0],DVopt[1],DVopt[2],DVopt[3],DVopt[4],DVopt[5],DVopt[6],DVopt[7])    

    # print("xf: ", repr(xf))
    # print("yf: ", repr(yf))
    # print("nturbs: ", len(xf))

