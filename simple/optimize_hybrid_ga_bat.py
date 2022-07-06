import numpy as np
from shapely.geometry import Polygon, Point
from gradient_free import GeneticAlgorithm, GreedyAlgorithm
from model_hybrid import SimpleHybrid
import floris.tools as wfct
import matplotlib.pyplot as plt
from hybrid_optimization_tools import hybrid_sweep, hybrid_sweep_w_battery
import time


def variables_to_chromosome(bits,bounds,variable_type,variable_array):  
        """convert the binary chromosomes to design variable values"""      

        nbits = sum(bits)
        nvars = len(bits)
        chromosome = np.zeros(nbits)

        start_bit = 0
        for i in range(nvars):
            if variable_type[i] == "int":
                low_val = bounds[i][0]
                int_val = variable_array[i]-low_val
                sum_val = 0
                for j in range(bits[i]-1,-1,-1):
                    if sum_val + 2**j <= int_val:
                        chromosome[start_bit+j] = 1
                        sum_val += 2**j
            
            elif variable_type[i] == "float":
                low_val = bounds[i][0]
                range_val = bounds[i][1]-bounds[i][0]
                float_val = variable_array[i]-low_val
                ndiscretizations = 2**bits[i]
                discretization = range_val/(ndiscretizations-1)
                sum_val = 0
                for j in range(bits[i]-1,-1,-1):
                    if (sum_val + 2**j)*discretization <= float_val:
                        chromosome[start_bit+j] = 1
                        sum_val += 2**j

            start_bit += bits[i]
        
        return chromosome


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


def create_plant_array(turbine_x,turbine_y,solar_x,solar_y):

    global xlocs
    global ylocs

    plant_array = np.zeros(len(xlocs),dtype=int)
    nturbs = len(turbine_x)
    narrays = len(solar_x)

    for i in range(len(xlocs)):
        for j in range(nturbs):
            if xlocs[i] == turbine_x[j] and ylocs[i] == turbine_y[j]:
                plant_array[i] = 1
        for j in range(narrays):
            if xlocs[i] == solar_x[j] and ylocs[i] == solar_y[j]:
                plant_array[i] = 2

    return plant_array


def capacity_obj_w_costs(inputs):

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

    global min_energy
    global solar_cost
    global battery_cost
    global under_penalty

    global function_calls

    global fixed_turbine_x
    global fixed_turbine_y
    global fixed_solar_x
    global fixed_solar_y

    # plant_array = inputs[0:-1]
    # battery_storage = inputs[-1]
    plant_array = inputs[:]
    battery_storage = 15.0

    t_x, t_y, s_x, s_y = create_farm(plant_array)
    turbine_x = np.append(fixed_turbine_x,t_x)
    turbine_y = np.append(fixed_turbine_y,t_y)
    solar_x = np.append(fixed_solar_x,s_x)
    solar_y = np.append(fixed_solar_y,s_y)

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
    # solar_solar = plant.get_solar_solar_distance()
    solar_matrix = plant.get_solar_solar_distance_matrix()
    solar_solar = np.ndarray.flatten(solar_matrix)

    if (len(solar_x)==0 and len(turbine_x)==0):
        return 1E20
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
        function_calls += 1
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

        wo_bat = np.zeros(nhours)  
        w_bat = np.zeros(nhours)
        for i in range(nhours):
            wo_bat[i] = sum(real_energy[0:i])
            w_bat[i] = sum(real_energy_with_battery[0:i])
        
        sc = sum(plant.solar_capacity_mw)
        wc = len(turbine_x)*2.5
        exponent = -0.00174
        # exponent = -0.0
        cost_solar = sc*solar_cost
        cost_wind = wc
        cost_battery = battery_cost*battery_storage

        if len(solar_x) < 2:
            solar_cable = 0.0
        else:
            solar_cable = 0.0
            for i in range(len(solar_x)):
                solar_cable += min(solar_matrix[i][:])
        # cost_cable = solar_cable/1000.0
        # total_cost = (cost_solar*0.9 + cost_cable*0.1 +cost_wind)*(2./3.+1./3.*np.exp(exponent*(wc+sc)**2))
        total_cost = (cost_solar + cost_wind)*(2./3.+1./3.*np.exp(exponent*(wc+sc)**2)) + penalty_cost + cost_battery
        # print(real_energy_with_battery)
        return total_cost/(np.sum(real_energy_with_battery))


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

    global fixed_turbine_x
    global fixed_turbine_y
    global fixed_solar_x
    global fixed_solar_y

    function_calls = 0

    battery_cost = 0.75
    under_penalty = 0.5

    solar_cost = 1.0
    min_energy = 10.0
    interconnect = 25.0

    np.random.seed(seed=1)

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

    side = 1600.0
    grid_size = 10
    x_grid = np.linspace(0.0,side,grid_size)
    y_grid = np.linspace(0.0,side,grid_size)

    solar_width = x_grid[1]-x_grid[0]-1.0
    solar_height = y_grid[1]-y_grid[0]-1.0

    xlocs, ylocs = np.meshgrid(x_grid,y_grid)
    xlocs = np.ndarray.flatten(xlocs)
    ylocs = np.ndarray.flatten(ylocs)
    
    plant = SimpleHybrid()
    plant.floris_model = wfct.floris_interface.FlorisInterface("/Users/astanley/Data/turbines/2_5mw_118d_88h.json")
    plant.floris_model.set_gch(False)

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

    battery_size = 0.0
    # gradient-free optimization
    start_time = time.time()

    fixed_tx0 = np.array([])
    fixed_tx1 = np.array([])
    fixed_tx2 = np.array([])
    fixed_tx3 = np.array([])
    fixed_ty0 = np.array([])
    fixed_ty1 = np.array([])
    fixed_ty2 = np.array([])
    fixed_ty3 = np.array([])
    fixed_sx0 = np.array([])
    fixed_sx1 = np.array([])
    fixed_sx2 = np.array([])
    fixed_sx3 = np.array([])
    fixed_sy0 = np.array([])
    fixed_sy1 = np.array([])
    fixed_sy2 = np.array([])
    fixed_sy3 = np.array([])

    nsections = 4
    converge_counter = 0
    total_iters = 0

    section = 0

    best_sol = 1E6

    while converge_counter < nsections:

        total_iters += 1

        if section == 0:
            xlocs, ylocs = np.meshgrid(x_grid[0:int(grid_size/2)],y_grid[0:int(grid_size/2)])
            fixed_turbine_x = np.concatenate((fixed_tx1,fixed_tx2,fixed_tx3))
            fixed_turbine_y = np.concatenate((fixed_ty1,fixed_ty2,fixed_ty3))
            fixed_solar_x = np.concatenate((fixed_sx1,fixed_sx2,fixed_sx3))
            fixed_solar_y = np.concatenate((fixed_sy1,fixed_sy2,fixed_sy3))

            xlocs = np.ndarray.flatten(xlocs)
            ylocs = np.ndarray.flatten(ylocs)

            start_array = create_plant_array(fixed_tx0,fixed_ty0,fixed_sx0,fixed_sy0)

        elif section == 1:
            xlocs, ylocs = np.meshgrid(x_grid[int(grid_size/2):],y_grid[0:int(grid_size/2)])
            fixed_turbine_x = np.concatenate((fixed_tx0,fixed_tx2,fixed_tx3))
            fixed_turbine_y = np.concatenate((fixed_ty0,fixed_ty2,fixed_ty3))
            fixed_solar_x = np.concatenate((fixed_sx0,fixed_sx2,fixed_sx3))
            fixed_solar_y = np.concatenate((fixed_sy0,fixed_sy2,fixed_sy3))

            xlocs = np.ndarray.flatten(xlocs)
            ylocs = np.ndarray.flatten(ylocs)

            start_array = create_plant_array(fixed_tx1,fixed_ty1,fixed_sx1,fixed_sy1)

        elif section == 2:
            xlocs, ylocs = np.meshgrid(x_grid[0:int(grid_size/2)],y_grid[int(grid_size/2):])
            fixed_turbine_x = np.concatenate((fixed_tx0,fixed_tx1,fixed_tx3))
            fixed_turbine_y = np.concatenate((fixed_ty0,fixed_ty1,fixed_ty3))
            fixed_solar_x = np.concatenate((fixed_sx0,fixed_sx1,fixed_sx3))
            fixed_solar_y = np.concatenate((fixed_sy0,fixed_sy1,fixed_sy3))

            xlocs = np.ndarray.flatten(xlocs)
            ylocs = np.ndarray.flatten(ylocs)

            start_array = create_plant_array(fixed_tx2,fixed_ty2,fixed_sx2,fixed_sy2)
            
        elif section == 3:
            xlocs, ylocs = np.meshgrid(x_grid[int(grid_size/2):],y_grid[int(grid_size/2):])
            fixed_turbine_x = np.concatenate((fixed_tx0,fixed_tx1,fixed_tx2))
            fixed_turbine_y = np.concatenate((fixed_ty0,fixed_ty1,fixed_ty2))
            fixed_solar_x = np.concatenate((fixed_sx0,fixed_sx1,fixed_sx2))
            fixed_solar_y = np.concatenate((fixed_sy0,fixed_sy1,fixed_sy2))

            xlocs = np.ndarray.flatten(xlocs)
            ylocs = np.ndarray.flatten(ylocs)

            start_array = create_plant_array(fixed_tx3,fixed_ty3,fixed_sx3,fixed_sy3)
        
        ga = GeneticAlgorithm()

        nvars = len(xlocs)
        ga.bits = np.zeros(nvars,dtype=int)
        ga.bounds = np.zeros((nvars+1,2))
        ga.variable_type = np.array([])
        for i in range(int(grid_size/2*grid_size/2)):
            ga.bits[i] = 2
            ga.variable_type = np.append(ga.variable_type,"int")
            ga.bounds[i] = (0,3)

        # ga.bits = np.append(ga.bits,8)
        # ga.variable_type = np.append(ga.variable_type,"float")
        # ga.bounds[-1,:] = (0.0,60.0)

        ga.objective_function = capacity_obj_w_costs

        ga.population_size = 25
        ga.convergence_iters = 10
        ga.max_generation = 1000

        ga.crossover_rate = 0.1
        ga.mutation_rate = 0.01

        # ga.optimize_ga(initialize="limit",crossover="chunk")
        # start_individual = np.append(start_array,battery_size)
        start_individual = np.copy(start_array)
        start_chromosome = variables_to_chromosome(ga.bits,ga.bounds,ga.variable_type,start_individual)
        print(start_chromosome)
        ga.optimize_ga(initialize="limit",start_individual=start_chromosome)

        opt_val = ga.optimized_function_value
        DVopt = ga.optimized_design_variables
        np.set_printoptions(precision=16)
        print(DVopt)
        if opt_val < best_sol:

            best_sol = opt_val
            t_x, t_y, s_x, s_y = create_farm(DVopt[0:-1])

            if section == 0:
                fixed_tx0 = np.copy(t_x)
                fixed_ty0 = np.copy(t_y)
                fixed_sx0 = np.copy(s_x)
                fixed_sy0 = np.copy(s_y)
            elif section == 1:
                fixed_tx1 = np.copy(t_x)
                fixed_ty1 = np.copy(t_y)
                fixed_sx1 = np.copy(s_x)
                fixed_sy1 = np.copy(s_y)
            if section == 2:
                fixed_tx2 = np.copy(t_x)
                fixed_ty2 = np.copy(t_y)
                fixed_sx2 = np.copy(s_x)
                fixed_sy2 = np.copy(s_y)
            if section == 3:
                fixed_tx3 = np.copy(t_x)
                fixed_ty3 = np.copy(t_y)
                fixed_sx3 = np.copy(s_x)
                fixed_sy3 = np.copy(s_y)
        
        else:
            converge_counter += 1

        tx = np.concatenate((fixed_tx0,fixed_tx1,fixed_tx2,fixed_tx3))
        ty = np.concatenate((fixed_ty0,fixed_ty1,fixed_ty2,fixed_ty3))
        sx = np.concatenate((fixed_sx0,fixed_sx1,fixed_sx2,fixed_sx3))
        sy = np.concatenate((fixed_sy0,fixed_sy1,fixed_sy2,fixed_sy3))
        plot_hybrid(tx,ty,sx,sy)
        print("section finished")
        # battery_size = DVopt[-1]
        print("best_sol: ", best_sol)
        print("tx: ", tx)
        print("ty: ", ty)
        print("sx: ", sx)
        print("sy: ", sy)
        # print("battery_size: ", battery_size)

        section = (section+1)%nsections
    

    run_time = time.time()-start_time

    save_filename = "ga_1.txt"
    file = open('%s'%save_filename, 'w')
    file.write("SOLUTION" + '\n' + '\n')
    file.write("optimal solution: " + '%s'%best_sol + '\n')
    file.write("optimal array: " + '%s'%DVopt + '\n')
    file.write("turbine_x = np." + '%s'%repr(t_x) + '\n')
    file.write("turbine_y = np." + '%s'%repr(t_y) + '\n')
    file.write("nturbs = " + '%s'%(len(t_x)) + '\n')
    file.write("solar_x = np." + '%s'%repr(s_x) + '\n')
    file.write("solar_y = np." + '%s'%repr(s_y) + '\n')
    file.write("nsolar = " + '%s'%(len(s_x)) + '\n')
    file.write("function_calls = " + '%s'%(function_calls) + '\n')
    # file.write("battery_storage = " + '%s'%(battery_storage) + '\n')
    file.write("run_time = " + '%s'%(run_time) + '\n' + '\n' + '\n' + '\n')

    file.write("SETUP" + '\n' + '\n')
    file.write("battery_cost = " + '%s'%(battery_cost) + '\n')
    file.write("solar_cost = " + '%s'%(solar_cost) + '\n')
    file.write("under_penalty = " + '%s'%(under_penalty) + '\n')
    file.write("min_energy = " + '%s'%(min_energy) + '\n')
    file.write("interconnect = " + '%s'%(interconnect) + '\n')

    file.write("side = " + '%s'%(side) + '\n')
    file.write("grid_size = " + '%s'%(grid_size) + '\n')

    file.write("plant.solar_capacity_factor = " + '%s'%(plant.solar_capacity_factor) + '\n')
    file.write("plant.acres_per_MW = " + '%s'%(plant.acres_per_MW) + '\n')
    file.write("plant.shadow_x = " + '%s'%(plant.shadow_x) + '\n')
    file.write("plant.shadow_y = " + '%s'%(plant.shadow_y) + '\n')
    file.write("plant.shadow_scale = " + '%s'%(plant.shadow_scale) + '\n')

    file.write("ws_array = np." + '%s'%repr(ws_array) + '\n')
    file.write("sr_array = np." + '%s'%repr(sr_array) + '\n')
    file.write("wd_array = np." + '%s'%repr(wd_array) + '\n')

    file.write("D = " + '%s'%(D) + '\n' + '\n')
    file.write("wind_wind_spacing = " + '%s'%(wind_wind_spacing) + '\n')
    file.write("solar_solar_spacing = " + '%s'%(solar_solar_spacing) + '\n')
    file.write("wind_solar_spacing = " + '%s'%(wind_solar_spacing) + '\n')

    file.close()

    # plot_hybrid(plant_array)
    plt.show()


