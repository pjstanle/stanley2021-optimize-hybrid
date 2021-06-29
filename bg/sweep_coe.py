import numpy as np
from shapely.geometry import Polygon, Point, LineString, MultiPolygon
from gradient_free import GeneticAlgorithm
from model_hybrid import SimpleHybrid
import floris.tools as wfct
import matplotlib.pyplot as plt
import time
from shapely.ops import unary_union
from boundary_grid import makeBoundary, makeGrid


def calculate_flicker(turbine_x,turbine_y,solar_poly):

        global rotor_diameter
        global plant

        def _gaus2d(sx,sy,dx,dy):
            return np.exp(-(dx**2/(2.0*sx**2)+dy**2/(2.0*sy**2)))

        nturbs = len(turbine_x)

        sigma_x = plant.shadow_x*rotor_diameter
        sigma_y = plant.shadow_y*rotor_diameter
        solar_loc = solar_poly.centroid
        shadow_multiplier = 1.0
        for j in range(nturbs):
            dx = turbine_x[j] - solar_loc.x
            dy = turbine_y[j] - solar_loc.y
            mult = 1.0 - _gaus2d(sigma_x,sigma_y,dx,dy)
            shadow_multiplier *= mult 
        
        return shadow_multiplier


def get_solar_polygon(solar_x,solar_y,solar_width,solar_height):

    solar_pts = [(solar_x-solar_width/2.0,solar_y-solar_height/2.0),
                (solar_x+solar_width/2.0,solar_y-solar_height/2.0),
                (solar_x+solar_width/2.0,solar_y+solar_height/2.0),
                (solar_x-solar_width/2.0,solar_y+solar_height/2.0),
                (solar_x-solar_width/2.0,solar_y-solar_height/2.0)]

    
    solar_poly = Polygon(solar_pts)
    return solar_poly


def combine_solar_polygons(poly0,poly1,poly2):
    nsolar = 0
    if poly0.area > 0.0:
        nsolar += 1
    if poly1.area > 0.0:
        nsolar += 1
    if poly2.area > 0.0:
        nsolar += 1
    
    solar_polys = [0]*nsolar
    counter = 0
    if poly0.area > 0.0:
        solar_polys[counter] = poly0
        counter += 1
    if poly1.area > 0.0:
        solar_polys[counter] = poly1
        counter += 1
    if poly2.area > 0.0:
        solar_polys[counter] = poly2
        counter += 1
    
    polygons = [Point(i, 0).buffer(0.7) for i in range(5)]
    solar_union = unary_union(polygons)
    solar_union = unary_union(solar_polys)
    if type(solar_union) == Polygon:
        solar_union = MultiPolygon([solar_union])

    return solar_union


def delete_solar_turbines(initial_turbine_x,initial_turbine_y,solar_union):
    # delete turbines close to solar
    keep_turbines = np.ones(len(initial_turbine_x),dtype=int)
    for i in range(len(initial_turbine_x)):
        for j in range(len(solar_union)):
            distance = get_wind_solar_distance(initial_turbine_x[i],initial_turbine_y[i],solar_union[j])
            if distance < wind_solar_spacing:
                keep_turbines[i] = 0
    
    turbine_x = np.zeros(np.sum(keep_turbines))
    turbine_y = np.zeros(np.sum(keep_turbines))
    counter = 0
    for i in range(len(initial_turbine_x)):
        if keep_turbines[i] == 1:
            turbine_x[counter] = initial_turbine_x[i]
            turbine_y[counter] = initial_turbine_y[i]
            counter += 1
    
    return turbine_x, turbine_y


def get_wind_solar_distance(turbine_x,turbine_y,solar_poly):
              
    global rotor_diameter
    
    x,y = solar_poly.boundary.coords.xy
    solar_pts = np.zeros((len(x),2))
    solar_pts[:,0] = x[:]
    solar_pts[:,1] = y[:]

    solar_line = LineString(solar_pts)
    turbine_point = Point(turbine_x,turbine_y)
    spacing = solar_line.distance(turbine_point)
    if solar_poly.contains(turbine_point):
        spacing *= -1.0

    return spacing-rotor_diameter


def plot_hybrid(turbine_x,turbine_y,solar_polys):

    global boundary_poly
    global rotor_diameter

    plt.cla()

    # plot boundaries
    bx,by = boundary_poly.boundary.coords.xy
    bx = np.append(bx,bx[0])
    by = np.append(by,by[0])
    plt.plot(bx,by,"--k",linewidth=0.5)

    # plot solar
    try:
        for j in range(len(solar_polys)):
            x,y = solar_polys[j].boundary.coords.xy
            plt.fill(x,y,color="C1")
    except:
        x,y = solar_polys.boundary.coords.xy
        plt.fill(x,y,color="C1")

    # plot turbines
    for i in range(len(turbine_x)):
        turb = plt.Circle((turbine_x[i],turbine_y[i]),rotor_diameter/2.0,facecolor="C0")
        plt.gca().add_patch(turb)

    
    plt.axis("equal")
    plt.pause(0.0001)


def obj_COE(inputs):

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

    global function_calls

    global boundary_poly
    global operation

    global best_answer

    # turbine variables
    boundary_start = inputs[0]
    boundary_turbs = inputs[1]
    grid_rows = inputs[2]
    grid_cols = inputs[3]
    grid_spacing_x = inputs[4]
    grid_spacing_y = inputs[5]
    grid_shear = inputs[6]
    grid_rotation = inputs[7]
    grid_center_x = inputs[8]
    grid_center_y = inputs[9]
    grid_boundary_setback = inputs[10]

    # solar variables
    solar_x = inputs[11:14]
    solar_y = inputs[14:17]
    solar_width = inputs[17:20]
    solar_height = inputs[20:23]

    # battery variables
    battery_storage = inputs[23]

    # assign turbine variables
    boundary_turbine_x, boundary_turbine_y = makeBoundary(boundary_start,boundary_turbs,boundary_poly)
    grid_turbine_x, grid_turbine_y = makeGrid(grid_rows,grid_cols,grid_spacing_x,grid_spacing_y,grid_shear,
                                                grid_rotation,grid_center_x,grid_center_y,grid_boundary_setback,boundary_poly)
    initial_turbine_x = np.append(boundary_turbine_x,grid_turbine_x)
    initial_turbine_y = np.append(boundary_turbine_y,grid_turbine_y)

    poly0 = get_solar_polygon(solar_x[0],solar_y[0],solar_width[0],solar_height[0]).intersection(boundary_poly)
    poly1 = get_solar_polygon(solar_x[1],solar_y[1],solar_width[1],solar_height[1]).intersection(boundary_poly)
    poly2 = get_solar_polygon(solar_x[2],solar_y[2],solar_width[2],solar_height[2]).intersection(boundary_poly)

    solar_union = combine_solar_polygons(poly0,poly1,poly2)

    turbine_x, turbine_y = delete_solar_turbines(initial_turbine_x,initial_turbine_y,solar_union)
    
    plant.turbine_x = turbine_x
    plant.turbine_y = turbine_y


    # check distances
    wind_wind = plant.get_wind_wind_distance()
    if len(wind_wind) > 0:
        min_wind_wind = min(wind_wind)
    else:
        min_wind_wind = 1E6


    if min_wind_wind < wind_wind_spacing:
        return 1E20
    elif len(solar_union) > 1:
        return 1E20
    else:
        function_calls += 1
        smpa = 4046.86 # square meters per acre
        solar_capacity = solar_union.area/smpa/plant.acres_per_MW
        wind_capacity = len(turbine_x)*2.5

        nhours = len(wd_array)
        gross_energy = np.zeros(nhours)
        
        shadow_scale = 0.5

        array_capacity = np.zeros(len(solar_union))
        shadow_multiplier = np.zeros(len(solar_union))
        for j in range(len(solar_union)):
            array_capacity[j] = solar_union[j].area/smpa/plant.acres_per_MW
            shadow_multiplier[j] = calculate_flicker(turbine_x,turbine_y,solar_union[j])

        for i in range(nhours):
            wind_energy = plant.get_wind_power(wd_array[i],ws_array[i])
            solar_energy = 0.0
            for j in range(len(solar_union)):
                solar_energy += array_capacity[j]*sr_array[i]*shadow_multiplier[j]**shadow_scale
            gross_energy[i] = wind_energy + solar_energy

        real_energy = np.zeros(nhours)
        battery_charge = np.zeros(nhours+1)
        real_energy_with_battery = np.zeros(nhours)
        
        battery_charge[0] = battery_storage

        if operation == "store_curtailment":
            for i in range(nhours):
                if gross_energy[i] > interconnect:
                    real_energy[i] = interconnect
                else:
                    real_energy[i] = gross_energy[i]
                curtailment = gross_energy[i]-real_energy[i]

                # charge the battery
                if curtailment > 0.0 and battery_charge[i] < battery_storage:
                    battery_charge[i+1] = battery_charge[i]+curtailment
                    real_energy_with_battery[i] = real_energy[i]
                    if battery_charge[i+1] > battery_storage:
                        battery_charge[i+1] = battery_storage
                # discharge the battery
                elif real_energy[i] < min_energy:
                    required_energy = min_energy-real_energy[i]
                    provided_energy = min(required_energy,battery_charge[i])
                    real_energy_with_battery[i] = real_energy[i]+provided_energy
                    battery_charge[i+1] = battery_charge[i]-provided_energy
                # do nothing with the battery
                else:
                    battery_charge[i+1] = battery_charge[i]
                    real_energy_with_battery[i] = real_energy[i]
            

        elif operation == "always_stay_charged":
            for i in range(nhours):
                # potentially charge the battery
                if gross_energy[i] > min_energy:
                    # if battery is already charged
                    if battery_charge[i] == battery_storage:
                        if gross_energy[i] > interconnect:
                            real_energy[i] = interconnect
                        else:
                            real_energy[i] = gross_energy[i]
                        real_energy_with_battery[i] = real_energy[i]
                        battery_charge[i+1] = battery_charge[i]
                    # if battery needs to be charged
                    else:
                        charge_needed = battery_storage-battery_charge[i]
                        charge_available = gross_energy[i]-min_energy
                        charge_provided = min(charge_needed,charge_available)
                        battery_charge[i+1] = battery_charge[i] + charge_provided
                        real_energy[i] = gross_energy[i]-charge_provided
                        if real_energy[i] > interconnect:
                            real_energy[i] = interconnect
                        real_energy_with_battery[i] = real_energy[i]
                        
                # potentially discharge the battery
                else:
                    real_energy[i] = gross_energy[i]
                    required_energy = min_energy-real_energy[i]
                    provided_energy = min(required_energy,battery_charge[i])
                    real_energy_with_battery[i] = real_energy[i]+provided_energy
                    battery_charge[i+1] = battery_charge[i]-provided_energy

        exponent = -0.002

        cost_solar = solar_capacity*solar_cost*(2./3.+1./3.*np.exp(exponent*(solar_capacity/2.5)**2))
        cost_wind = wind_capacity
        cost_battery = battery_storage*battery_cost

        # cost_solar = cost_solar*(2./3.+1./3.*np.exp(exponent*(solar_capacity/2.5)**2))
        cost_wind = cost_wind*(2./3.+1./3.*np.exp(exponent*(wind_capacity/2.5)**2))
        cost_battery = cost_battery*(2./3.+1./3.*np.exp(exponent*(battery_storage/2.5)**2))

        total_cost = cost_solar + cost_wind + cost_battery
        total_energy = (np.sum(real_energy_with_battery))

        coe_obj = total_cost/total_energy
        
        return coe_obj, total_energy, (2./3.+1./3.*np.exp(exponent*(solar_capacity/2.5)**2)), solar_capacity,solar_union.area



if __name__=="__main__":

    global rotor_diameter
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
    global operation
    global boundary_poly
    global best_answer

    best_answer = 1E20

    operation = "always_stay_charged"
    # operation = "store_curtailment"
    function_calls = 0

    battery_cost = 0.4
    solar_cost = 0.7
    min_energy = 0.0
    interconnect = 20.0
    under_penalty = 0.5

    sc = np.linspace(0.1,1.5,10)

    boundary_poly = Polygon(([0,0],[1000,-100],[1200,1000],[500,1100],[-200,1000])).buffer(200)

    # np.random.seed(seed=1)

    ws_array = np.array([0.95,1.00,1.00,0.95,0.85,0.7,0.5,0.3,0.15,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.15,0.3,0.5,0.7,0.85])*6.0+6.0
    sr_array = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.15,0.3,0.5,0.7,0.85,0.95,1.00,1.00,0.95,0.85,0.7,0.5,0.3,0.15,0.05,0.0,0.0])
    ws_array = np.append(ws_array,ws_array)
    sr_array = np.append(sr_array,sr_array)

    ws_array = ws_array - np.random.rand(len(ws_array))
    sr_array = sr_array - np.random.rand(len(sr_array))*0.1
    for i in range(len(sr_array)):
        if sr_array[i] < 0.0:
            sr_array[i] = 0.0

    wd_array = np.ones(len(ws_array))*270.0 + np.random.rand(len(ws_array))*30.0
    
    plant = SimpleHybrid()
    plant.floris_model = wfct.floris_interface.FlorisInterface("/Users/astanley/Data/turbines/2_5mw_118d_88h.json")
    rotor_diameter = plant.floris_model.floris.farm.turbines[0].rotor_diameter
    plant.acres_per_MW = 2.5
    plant.shadow_x = 2.0 # rotor diameters of the sigma of shadow loss in x (east west)
    plant.shadow_y = 0.5 # rotor diameters of the sigma of shadow loss in y (north south)

    D = 118.0
    wind_wind_spacing = D*4.0
    wind_solar_spacing = D*0.5


    inputs = np.zeros(24)
    
    inputs[11] = 500.0
    inputs[14] = 500.0
    inputs[23] = 0.0


    for k in range(10):
        # solar_cost = sc[k]
        solar_cost = 0.8
        interconnect = k*5

        power = np.zeros(100)
        coe = np.zeros(100)
        cost_thing = np.zeros(100)
        cap = np.zeros(100)
        area = np.zeros(100)

        size = np.linspace(0.0,2000.0,100)
        for i in range(100):
            inputs[17] = size[i]
            inputs[20] = size[i]
            coe[i],power[i],cost_thing[i],cap[i],area[i] = obj_COE(inputs)

        # plt.figure(2)
        # plt.title("aep")
        # plt.plot(area,power)
        
        # plt.figure(3)
        # plt.title("cost")
        # plt.plot(area,cost_thing)

        # plt.figure(4)
        # plt.title("cap")
        # plt.plot(area,cap)

        # plt.figure(5)
        # plt.title("area")
        # plt.plot(size,area)

        plt.figure(1)
        plt.title("coe")
        plt.plot(area,coe)
        plt.ylim(0,0.2)

    plt.show()