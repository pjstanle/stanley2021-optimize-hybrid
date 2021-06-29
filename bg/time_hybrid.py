import numpy as np
from shapely.geometry import Polygon, Point, LineString, MultiPolygon
from gradient_free import GeneticAlgorithm
from model_hybrid import SimpleHybrid
import floris.tools as wfct
import matplotlib.pyplot as plt
import time
from shapely.ops import unary_union
from boundary_grid import makeBoundary, makeGrid
import scipy.interpolate
from analyze_wind import init_wind_plant
from simple_shadow import SimpleShadow



def get_solar_polygon(solar_x,solar_y,solar_width,solar_height):

    solar_pts = [(solar_x-solar_width/2.0,solar_y-solar_height/2.0),
                (solar_x+solar_width/2.0,solar_y-solar_height/2.0),
                (solar_x+solar_width/2.0,solar_y+solar_height/2.0),
                (solar_x-solar_width/2.0,solar_y+solar_height/2.0),
                (solar_x-solar_width/2.0,solar_y-solar_height/2.0)]
    
    pts_x = np.array([solar_x-solar_width/2.0,solar_x+solar_width/2.0,solar_x+solar_width/2.0,solar_x-solar_width/2.0])
    pts_y = np.array([solar_y-solar_height/2.0,solar_y-solar_height/2.0,solar_y+solar_height/2.0,solar_y+solar_height/2.0])
    
    solar_poly = Polygon(solar_pts)
    return solar_poly, pts_x, pts_y


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


def objective_function(inputs):

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

    global wind_capex
    global solar_capex
    global battery_capex
    global HOPP_plant
    global ppa
    global turbine_rating
    global penalty_multiplier

    global time_array

    global shadow

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

    solar_x = inputs[11]
    solar_y = inputs[12]
    solar_width = inputs[13]
    solar_height = inputs[14]

    # battery variables
    battery_storage = inputs[15]

    s = time.time()
    # assign turbine variables
    boundary_turbine_x, boundary_turbine_y = makeBoundary(boundary_start,boundary_turbs,boundary_poly)
    grid_turbine_x, grid_turbine_y = makeGrid(grid_rows,grid_cols,grid_spacing_x,grid_spacing_y,grid_shear,
                                                grid_rotation,grid_center_x,grid_center_y,grid_boundary_setback,boundary_poly)
    initial_turbine_x = np.append(boundary_turbine_x,grid_turbine_x)
    initial_turbine_y = np.append(boundary_turbine_y,grid_turbine_y)

    poly,pts_x,pts_y = get_solar_polygon(solar_x,solar_y,solar_width,solar_height)
    poly = poly.intersection(boundary_poly)

   
    if type(poly)==Polygon:
        solar_union = MultiPolygon([poly])
    else:
        solar_union = poly

    turbine_x, turbine_y = delete_solar_turbines(initial_turbine_x,initial_turbine_y,solar_union)
    
    plant.turbine_x = turbine_x
    plant.turbine_y = turbine_y

    shadow.solar_x = pts_x
    shadow.solar_y = pts_y
    shadow.turbine_x = turbine_x
    shadow.turbine_y = turbine_y

    print("place turbines and solar panels: ", time.time()-s)

    
    s = time.time()
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
        print("check spacings: ", time.time()-s)

        s = time.time()
        function_calls += 1
        smpa = 4046.86 # square meters per acre
        solar_capacity = solar_union.area/smpa/plant.acres_per_MW
        wind_capacity = len(turbine_x)*turbine_rating

        if solar_capacity + wind_capacity > interconnect:
            return 1E6

        nhours = len(sr_array)
        gross_energy = np.zeros(nhours)
        
        array_capacity = np.zeros(len(solar_union))
        for j in range(len(solar_union)):
            array_capacity[j] = solar_union[j].area/smpa/plant.acres_per_MW
        print("calculate solar capacity: ", time.time()-s)

        s = time.time()
        if len(turbine_x) > 0:
            HOPP_plant.modify_coordinates(turbine_x,turbine_y)
            HOPP_plant.simulate(1)
            wind_energy = HOPP_plant.system_model.Outputs.gen
        else:
            wind_energy = np.zeros(len(sr_array))
        
        print("calculate wind energy: ", time.time()-s)

        shadow_time = 0.0
        solar_time = 0.0

        times = np.linspace(8,16,5)
        shadows = np.zeros(len(times))
        for i in range(len(times)):
            shadows[i] = shadow.calculate_losses(times[i], 180.0)

        shadow_loss_function = scipy.interpolate.interp1d(times, shadows, kind='cubic')

        for i in range(nhours):
            solar_energy = 0.0
            hour = time_array[i]%24
            # day = time_array[i]/24.0

            s = time.time()
            if hour > 8.0 and hour < 16.0 and solar_capacity > 0.0 and len(turbine_x) > 0:
                # shadow_loss = shadow.calculate_losses(hour, day)
                shadow_loss = shadow_loss_function(hour)
            else: 
                shadow_loss = 0.0
            shadow_time += (time.time()-s)

            s = time.time()
            for j in range(len(solar_union)):
                solar_energy += array_capacity[j]*sr_array[i]*(1.0-shadow_loss)
            solar_time += (time.time()-s)
            gross_energy[i] = wind_energy[i]/1000.0 + solar_energy

        print("shadow time: ", shadow_time)
        print("solar time: ", solar_time)
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

        cost_solar = solar_capex(solar_capacity)
        cost_wind = wind_capex(wind_capacity)
        cost_battery = battery_capex(battery_storage)
        cost_om = om_function(solar_capacity+wind_capacity+battery_storage)

        total_cost = cost_solar + cost_wind + cost_battery

        penalty = 0.0
        for i in range(len(real_energy_with_battery)):
            if real_energy_with_battery[i] < min_energy:
                penalty += (min_energy-real_energy_with_battery[i])*ppa[i]*penalty_multiplier

        # if penalty > 0.0:
        #     return 1E6

        fcr = 0.063
        annual_cost = fcr*total_cost + cost_om + penalty

        
        total_energy = np.sum(real_energy_with_battery)

        annual_income = np.sum(real_energy_with_battery*ppa)


        coe = annual_cost/total_energy
        aep = -total_energy
        pro = -(annual_income - annual_cost)/1E6
        
        obj = aep

        if obj < best_answer:
            best_answer = obj
            plot_hybrid(turbine_x,turbine_y,solar_union)

            print("cost solar: ", cost_solar/1E6)
            print("cost wind: ", cost_wind/1E6)
            print("cost_battery: ", cost_battery/1E6)
            print("cost_om: ", cost_om/1E6)
            print("capacity solar: ", solar_capacity)
            print("capacity wind: ", wind_capacity)
            print("capacity battery: ", battery_storage)
            print("annual cost: ", annual_cost/1E6)
            print("total energy: ", total_energy/1E3)
            print("penalty: ", penalty/1E6)
            print("coe: ", coe)
            print("aep: ", aep/1E6)
            print("pro: ", pro)

        return obj



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
    global penalty_multiplier
    global function_calls
    global operation
    global boundary_poly
    global best_answer
    global wind_capex
    global solar_capex
    global battery_capex
    global om_function
    global HOPP_plant
    global ppa
    global turbine_rating
    global time_array
    global shadow

    powercurve_filename = 'high_7r_200d_135h.txt'
    rotor_diameter = 200.0
    hub_height = 135.0
    turbine_rating = 7.0

    shadow = SimpleShadow()
    shadow.hub_height = hub_height
    shadow.rotor_diameter = rotor_diameter
    shadow.tower_width = 6.0
    shadow.latitude = 40.966

    HOPP_plant = init_wind_plant(hub_height,rotor_diameter,powercurve_filename)


    capex_cost = np.array([5*1786.0,2*1786.0,1786.0,1622.0,1528.0,1494.0,1470.0,1421.0,1408.0,1400.0]) # $/kW realistic
    capex_size = np.array([0.0,1.0,20.0,50.0,100.0,150.0,200.0,400.0,1000.0,100000.0]) # MW
    cost = capex_size*capex_cost*1000.0
    wind_capex = scipy.interpolate.interp1d(capex_size, cost, kind='cubic')
    solar_capex = scipy.interpolate.interp1d(capex_size, cost*0.75, kind='cubic')
    battery_capex = scipy.interpolate.interp1d(capex_size, cost*0.5, kind='cubic')

    def om_function(capacity):
        return 37.0*capacity*1000.0

    best_answer = 1E20

    operation = "always_stay_charged"
    function_calls = 0

    battery_cost = 0.4
    solar_cost = 1.0
    min_energy = 600.0
    interconnect = 600.0
    penalty_multiplier = 0.0

    boundary_poly = Polygon(([0,0],[1000,-1000],[1200,1000],[500,1100],[-200,1000])).buffer(2000)

    solar_array = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.15,0.3,0.5,0.7,0.85,0.95,1.00,1.00,0.95,0.85,0.7,0.5,0.3,0.15,0.05,0.0,0.0])
    sr_array = np.array([])
    for i in range(365):
        sr_array = np.append(sr_array,solar_array)
    
    np.random.seed(seed=1423)
    sr_array = sr_array - np.random.rand(len(sr_array))*0.1
    for i in range(len(sr_array)):
        if sr_array[i] < 0.0:
            sr_array[i] = 0.0

    data = HOPP_plant.site.wind_resource.data["data"]
    nhours = np.shape(data)[0]
    time_array = np.arange(nhours)

    def PPA_curve(hour):
        time = hour%24
        p1 = np.cos((8-time)/np.pi)*65.0
        p2 = np.cos((18-time)/np.pi)*65.0
        p3 = 55.0
        return max([p1,p2,p3])

    ppa = np.zeros(nhours)
    for i in range(nhours):
        ppa[i] = PPA_curve(time_array[i])

    
    plant = SimpleHybrid()
    plant.floris_model = wfct.floris_interface.FlorisInterface("/Users/astanley/Data/turbines/2_5mw_118d_88h_1jensen.json")
    rotor_diameter = plant.floris_model.floris.farm.turbines[0].rotor_diameter
    plant.acres_per_MW = 4.0
    plant.shadow_x = 2.0 # rotor diameters of the sigma of shadow loss in x (east west)
    plant.shadow_y = 0.5 # rotor diameters of the sigma of shadow loss in y (north south)

    D = rotor_diameter
    wind_wind_spacing = D*4.0
    wind_solar_spacing = D*0.5

    plant = SimpleHybrid()
    plant.acres_per_MW = 5.0
    plant.shadow_x = 2.0 # rotor diameters of the sigma of shadow loss in x (east west)
    plant.shadow_y = 0.5 # rotor diameters of the sigma of shadow loss in y (north south)

    D = rotor_diameter
    wind_wind_spacing = D*4.0
    wind_solar_spacing = D*0.5

    # gradient-free optimization
    start_time = time.time()

    # turbine variables
    boundary_start = 0.0
    boundary_turbs = 10
    grid_rows = 4
    grid_cols = 4
    grid_spacing_x = D*5
    grid_spacing_y = D*5
    grid_shear = 0.2
    grid_rotation = 0.4
    grid_center_x = 500
    grid_center_y = 500
    grid_boundary_setback = D*5

    solar_x = 500
    solar_y = 500
    solar_width = 1000
    solar_height = 400

    # battery variables
    battery_storage = 50

    inputs = np.zeros(16)
    inputs[0] = boundary_start
    inputs[1] = boundary_turbs
    inputs[2] = grid_rows
    inputs[3] = grid_cols
    inputs[4] = grid_spacing_x
    inputs[5] = grid_spacing_y
    inputs[6] = grid_shear
    inputs[7] = grid_rotation
    inputs[8] = grid_center_x
    inputs[9] = grid_center_y
    inputs[10] = grid_boundary_setback
    inputs[11] = solar_x
    inputs[12] = solar_y
    inputs[13] = solar_width
    inputs[14] = solar_height
    inputs[15] = battery_storage

    start = time.time()
    print("start")
    obj = objective_function(inputs)
    print("whole run: ", time.time()-start)
    print(obj)

    plt.show()

    

   