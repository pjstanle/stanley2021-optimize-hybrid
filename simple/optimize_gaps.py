import numpy as np
from shapely.geometry import Polygon, Point, MultiPolygon
from gradient_free import GeneticAlgorithm, GreedyAlgorithm
from model_hybrid import SimpleHybrid
import floris.tools as wfct
import matplotlib.pyplot as plt
from hybrid_optimization_tools import hybrid_sweep, hybrid_sweep_w_battery
import time
from shapely.ops import unary_union


def create_farm(plant_array):

    global xlocs
    global ylocs

    nturbs = np.count_nonzero(plant_array==1)
    turbine_x = np.zeros(nturbs)
    turbine_y = np.zeros(nturbs)

    turbine_ind = 0
    for i in range(len(plant_array)):
        if int(plant_array[i]) == 1:
            turbine_x[turbine_ind] = xlocs[i]
            turbine_y[turbine_ind] = ylocs[i]
            turbine_ind += 1

    return turbine_x, turbine_y


def create_solar_ellipse(turbine_x,turbine_y):
    global ex
    global ey
    global NPTS

    poly_array = np.zeros((NPTS,2))
    R = np.linspace(0,2*np.pi-2*np.pi/NPTS,NPTS)
    for i in range(NPTS):
        poly_array[i,0] = ex*np.cos(R[i]) + turbine_x
        poly_array[i,1] = ey*np.sin(R[i]) + turbine_y
    
    return Polygon(poly_array)


def create_solar_fill(turbine_x,turbine_y):
    global boundary_poly
    
    solar_poly = boundary_poly
    nturbs = len(turbine_x)

    for i in range(nturbs):
        turbine_ellipse = create_solar_ellipse(turbine_x[i],turbine_y[i])
        solar_poly = solar_poly.difference(turbine_ellipse)

    return solar_poly


def capacity_obj_w_costs(inputs):

    global plant
    global wind_wind_spacing
    global wind_solar_spacing

    global wd_array
    global ws_array
    global sr_array

    global interconnect

    global min_energy
    global solar_cost
    global battery_cost
    global under_penalty

    global function_calls
    global ex
    global ey
    global best_solution

    # start = time.time()
    plant_array = inputs[:]
    battery_storage = 25.0

    turbine_x, turbine_y = create_farm(plant_array)

    plant.turbine_x = turbine_x
    plant.turbine_y = turbine_y

    wind_wind = plant.get_wind_wind_distance()

    if len(wind_wind) > 0:
        if  min(wind_wind) < wind_wind_spacing:
            return 1E20
        
    solar_poly = create_solar_fill(turbine_x,turbine_y)
    solar_area = solar_poly.area
    smpa = 4046.86 # square meters per acre
    solar_capacity = solar_area/smpa/plant.acres_per_MW
    plant_capacity = solar_capacity + len(turbine_x)*2.5

    if plant_capacity > 10.0*interconnect:
        return 1E20
    else:
        function_calls += 1
        nhours = len(wd_array)
        gross_energy = np.zeros(nhours)
        real_energy = np.zeros(nhours)

        # start_plant = time.time()
        for i in range(nhours):
            # gross_energy[i] = plant.get_plant_power(wd_array[i],ws_array[i],sr_array[i])
            wind_energy = plant.get_wind_power(wd_array[i],ws_array[i])
            solar_energy = solar_capacity*sr_array[i]
            gross_energy[i] = wind_energy+solar_energy

            if gross_energy[i] > interconnect:
                real_energy[i] = interconnect
            else:
                real_energy[i] = gross_energy[i]
        # print("get plant power: ", time.time()-start_plant)
        # print("get plant power: ", time.time()-start)
        # start = time.time()

        curtailment = gross_energy-real_energy

        battery_charge = np.zeros(nhours+1)
        real_energy_with_battery = np.zeros(nhours)

        penalty_cost = 0.0
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

            if real_energy_with_battery[i] < min_energy:
                penalty_cost += under_penalty*(min_energy-real_energy_with_battery[i])
                # penalty_cost += under_penalty

        # print("battery loop: ", time.time()-start)
        # start = time.time()

        # wo_bat = np.zeros(nhours)  
        # w_bat = np.zeros(nhours)
        # for i in range(nhours):
        #     wo_bat[i] = sum(real_energy[0:i])
        #     w_bat[i] = sum(real_energy_with_battery[0:i])
        
        # print("cumulative sum loop: ", time.time()-start)
        # start = time.time()

        sc = solar_capacity
        wc = len(turbine_x)*2.5
        exponent = -0.00174
        # exponent = -0.0

        # n_solar_groups, solar_areas = union_solar(solar_x,solar_y)

        if solar_poly.type == 'Polygon':
            n_solar_groups = 1
        else:
            n_solar_groups = len(solar_poly)

        # print("solar union: ", time.time()-start)
        # start = time.time()

        cost_solar = sc*solar_cost
        cost_solar += (n_solar_groups-1)*0.025
        cost_wind = wc
        cost_battery = battery_cost*battery_storage

        # if len(solar_x) < 2:
        #     solar_cable = 0.0
        # else:
        #     solar_cable = 0.0
        #     for i in range(len(solar_x)):
        #         solar_cable += min(solar_matrix[i][:])
        # cost_cable = solar_cable/1000.0
        # total_cost = (cost_solar + cost_wind)*(2./3.+1./3.*np.exp(exponent*(wc+sc)**2)) + penalty_cost + cost_battery
        total_cost = cost_solar*(2./3.+1./3.*np.exp(exponent*(sc)**2)) + cost_wind*(2./3.+1./3.*np.exp(exponent*(wc)**2)) + penalty_cost + cost_battery
        
        # print("wrap up: ", time.time()-start)

        obj = total_cost/(np.sum(real_energy_with_battery))
        
        if obj < best_solution:
            print(real_energy_with_battery)
            best_solution = obj
            plot_hybrid(turbine_x,turbine_y)

        return obj


def get_xy(A):
    x = np.zeros(len(A))
    y = np.zeros(len(A))
    for i in range(len(A)):
        x[i] = A[i][0]
        y[i] = A[i][1]
    return x,y


def plot_solar(geom,ax):
    if geom.type == 'Polygon':
        exterior_coords = geom.exterior.coords[:]
        x,y = get_xy(exterior_coords)
        ax.fill(x,y,"C1")

        for interior in geom.interiors:
            interior_coords = interior.coords[:]
            x,y = get_xy(interior_coords)
            ax.fill(x,y,"white")

    elif geom.type == 'MultiPolygon':

        for part in geom:
            exterior_coords = part.exterior.coords[:]
            x,y = get_xy(exterior_coords)
            ax.fill(x,y,"C1")

            for interior in part.interiors:
                interior_coords = interior.coords[:]
                x,y = get_xy(interior_coords)
                ax.fill(x,y,"white")


def plot_hybrid(turbine_x,turbine_y):

    plt.cla()


    solar_poly = create_solar_fill(turbine_x,turbine_y)
    plot_solar(solar_poly,plt.gca())

    for i in range(len(turbine_x)):
        turb = plt.Circle((turbine_x[i],turbine_y[i]),59.0,facecolor="C0")
        plt.gca().add_patch(turb)

    plt.axis("equal")
    plt.pause(0.0001)


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

    global solar_cost
    global under_penalty
    global battery_cost

    global function_calls

    global ex
    global ey
    global NPTS
    global boundary_poly

    global best_solution

    best_solution = 1E15

    function_calls = 0

    battery_cost = 0.75
    under_penalty = 0.5

    solar_cost = 1.0
    min_energy = 25.0
    interconnect = 50.0

    np.random.seed(seed=10)

    ws_array = np.array([0.95,1.00,1.00,0.95,0.85,0.7,0.5,0.3,0.15,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.15,0.3,0.5,0.7,0.85])*6.0+4.0
    sr_array = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.15,0.3,0.5,0.7,0.85,0.95,1.00,1.00,0.95,0.85,0.7,0.5,0.3,0.15,0.05,0.0,0.0])
    ws_array = np.append(ws_array,ws_array)
    sr_array = np.append(sr_array,sr_array)

    # ws_array = np.array([1.00,0.95,0.7,0.3,0.15,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.15,0.3,0.7])*6.0+4.0
    # sr_array = np.array([0.0,0.0,0.0,0.15,0.3,0.5,0.7,0.85,0.95,1.00,1.00,0.95,0.85,0.7,0.5,0.3,0.15,0.0])

    ws_array = ws_array - np.random.rand(len(ws_array))
    sr_array = sr_array - np.random.rand(len(sr_array))*0.1
    for i in range(len(sr_array)):
        if sr_array[i] < 0.0:
            sr_array[i] = 0.0

    wd_array = np.ones(len(ws_array))*270.0 + np.random.rand(len(ws_array))*30.0

    grid_size = 10
    side = 118.0*(grid_size-1)
    x_grid = np.linspace(0.0,side,grid_size)
    y_grid = np.linspace(0.0,side,grid_size)

    solar_width = x_grid[1]-x_grid[0]-1.0
    solar_height = y_grid[1]-y_grid[0]-1.0

    xlocs, ylocs = np.meshgrid(x_grid,y_grid)
    xlocs = np.ndarray.flatten(xlocs)
    ylocs = np.ndarray.flatten(ylocs)
    
    plant = SimpleHybrid()
    plant.floris_model = wfct.floris_interface.FlorisInterface("/Users/astanley/Data/turbines/2_5mw_118d_88h_1jensen.json")

    plant.solar_capacity_factor = 0.2
    plant.acres_per_MW = 5.0
    plant.shadow_x = 2.0 # rotor diameters of the sigma of shadow loss in x (east west)
    plant.shadow_y = 0.5 # rotor diameters of the sigma of shadow loss in y (north south)
    plant.shadow_scale = 0.5

    D = 118.0
    wind_wind_spacing = D*4.0
    solar_solar_spacing = 0.0
    wind_solar_spacing = D*0.5

    nlocs = len(xlocs)
    ntech = 2

    # gradient-free optimization
    start_time = time.time()

    xlocs, ylocs = np.meshgrid(x_grid,y_grid)
    xlocs = np.ndarray.flatten(xlocs)
    ylocs = np.ndarray.flatten(ylocs)

    farm_width = max(xlocs)-min(xlocs)
    farm_height = max(ylocs)-min(ylocs)
    dx = x_grid[1]-x_grid[0]
    dy = y_grid[1]-y_grid[0]


    ex = 300
    ey = 100
    NPTS = 100

    boundary_poly = Polygon(([min(xlocs),min(ylocs)],[max(xlocs),min(ylocs)],[max(xlocs),max(ylocs)],[min(xlocs),max(ylocs)]))

    tx = np.array([0.0,100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0])
    ty = np.ones(len(tx))*500.0
    sp = create_solar_fill(tx,ty)

    ga = GeneticAlgorithm()

    nvars = len(xlocs)
    ga.bits = np.zeros(nvars,dtype=int)
    ga.bounds = np.zeros((nvars+1,2))
    ga.variable_type = np.array([])
    for i in range(int(grid_size*grid_size)):
        ga.bits[i] = 1
        ga.variable_type = np.append(ga.variable_type,"int")
        ga.bounds[i] = (0,2)

    ga.objective_function = capacity_obj_w_costs

    ga.population_size = 25
    ga.convergence_iters = 10
    ga.max_generation = 1000

    ga.crossover_rate = 0.1
    ga.mutation_rate = 0.01

    ga.optimize_ga(initialize="limit")

    opt_val = ga.optimized_function_value
    DVopt = ga.optimized_design_variables
    np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})
    np.set_printoptions(precision=16)
    print(repr(DVopt))

    print("finished")
    print("best_sol: ", opt_val)
    # print("tx: ", tx)
    # print("ty: ", ty)
    # print("sx: ", sx)
    # print("sy: ", sy)
    

    run_time = time.time()-start_time
    plt.show()

   