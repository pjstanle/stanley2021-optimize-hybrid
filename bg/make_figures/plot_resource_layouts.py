
import numpy as np
from shapely.geometry import Polygon, Point, LineString, MultiPolygon
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


def get_xy(A):
    x = np.zeros(len(A))
    y = np.zeros(len(A))
    for i in range(len(A)):
        x[i] = A[i][0]
        y[i] = A[i][1]
    return x,y


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
    # try:
    #     for j in range(len(solar_polys)):
    #         x,y = solar_polys[j].boundary.coords.xy
    #         plt.fill(x,y,color="C1")
    # except:
    #     x,y = solar_polys.boundary.coords.xy
    #     plt.fill(x,y,color="C1")
    if type(solar_polys) == Polygon:
        x,y = solar_polys.boundary.coords.xy
        plt.fill(x,y,color="C1")
    else:
        for part in solar_polys:
            exterior_coords = part.exterior.coords[:]
            x,y = get_xy(exterior_coords)
            plt.fill(x,y,color="C1")

            for interior in part.interiors:
                interior_coords = interior.coords[:]
                x,y = get_xy(interior_coords)
                plt.fill(x,y,"white")


    # plot turbines
    for i in range(len(turbine_x)):
        # turb = plt.Circle((turbine_x[i],turbine_y[i]),rotor_diameter/2.0,facecolor="C0")
        turb = plt.Circle((turbine_x[i],turbine_y[i]),rotor_diameter/2.0,facecolor="C0")
        plt.gca().add_patch(turb)



def get_dvs(inputs):

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

    return turbine_x, turbine_y, solar_union


def objective_function(inputs):

    global plant
    global wind_wind_spacing
    global wind_solar_spacing

    global sr_array

    global interconnect
    global absolute_min_energy
    global penalty_min_energy

    global boundary_poly
    global operation

    global wind_capex
    global solar_capex
    global battery_capex
    global HOPP_plant
    global ppa
    global turbine_rating
    global penalty_multiplier

    global time_array

    global shadow

    global max_capacity

    global objective


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

    smpa = 4046.86 # square meters per acre
    solar_capacity = solar_union.area/smpa/plant.acres_per_MW
    wind_capacity = len(turbine_x)*turbine_rating

    if solar_capacity + wind_capacity > max_capacity:
        return 1E6

    nhours = len(sr_array)
    gross_energy = np.zeros(nhours)
    
    array_capacity = np.zeros(len(solar_union))
    for j in range(len(solar_union)):
        array_capacity[j] = solar_union[j].area/smpa/plant.acres_per_MW

    if len(turbine_x) > 0:
        HOPP_plant.modify_coordinates(turbine_x,turbine_y)
        HOPP_plant.simulate(1)
        wind_energy = HOPP_plant.system_model.Outputs.gen
    else:
        wind_energy = np.zeros(len(sr_array))

    times = np.linspace(8,16,5)
    shadows = np.zeros(len(times))
    if solar_capacity > 0.0 and len(turbine_x) > 0:
        for i in range(len(times)):
            shadows[i] = shadow.calculate_losses(times[i], 180.0)

    shadow_loss_function = scipy.interpolate.interp1d(times, shadows, kind='cubic')

    for i in range(nhours):
        solar_energy = 0.0
        hour = time_array[i]%24

        s = time.time()
        if hour > 8.0 and hour < 16.0:
            # shadow_loss = shadow.calculate_losses(hour, day)
            shadow_loss = shadow_loss_function(hour)
        else: 
            shadow_loss = 0.0

        for j in range(len(solar_union)):
            solar_energy += array_capacity[j]*sr_array[i]*(1.0-shadow_loss)
        gross_energy[i] = wind_energy[i]/1000.0 + solar_energy

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
            elif real_energy[i] < penalty_min_energy:
                required_energy = penalty_min_energy-real_energy[i]
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
            if gross_energy[i] > interconnect:
                # if battery is already charged
                if battery_charge[i] == battery_storage:
                    real_energy[i] = interconnect
                    real_energy_with_battery[i] = real_energy[i]
                    battery_charge[i+1] = battery_charge[i]
                # if battery needs to be charged
                else:
                    charge_needed = battery_storage-battery_charge[i]
                    charge_available = gross_energy[i]-interconnect
                    charge_provided = min(charge_needed,charge_available)
                    battery_charge[i+1] = battery_charge[i] + charge_provided
                    real_energy[i] = gross_energy[i]-charge_provided
                    if real_energy[i] > interconnect:
                        real_energy[i] = interconnect
                    real_energy_with_battery[i] = real_energy[i]

            elif gross_energy[i] > penalty_min_energy:
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
                    charge_available = gross_energy[i]-penalty_min_energy
                    charge_provided = min(charge_needed,charge_available)
                    battery_charge[i+1] = battery_charge[i] + charge_provided
                    real_energy[i] = gross_energy[i]-charge_provided
                    if real_energy[i] > interconnect:
                        real_energy[i] = interconnect
                    real_energy_with_battery[i] = real_energy[i]
                    
            # potentially discharge the battery
            else:
                real_energy[i] = gross_energy[i]
                required_energy = penalty_min_energy-real_energy[i]
                provided_energy = min(required_energy,battery_charge[i])
                real_energy_with_battery[i] = real_energy[i]+provided_energy
                battery_charge[i+1] = battery_charge[i]-provided_energy

    # if real_energy_with_battery.any() < absolute_min_energy:
    if min(real_energy_with_battery) < absolute_min_energy:
        return 1E6

    cost_solar = solar_capex(solar_capacity)
    cost_wind = wind_capex(wind_capacity)
    cost_battery = battery_capex(battery_storage)
    cost_om = om_function(solar_capacity+wind_capacity+battery_storage)

    total_cost = cost_solar + cost_wind + cost_battery

    penalty = 0.0
    if penalty_multiplier != 0.0:
        for i in range(len(real_energy_with_battery)):
            if real_energy_with_battery[i] < penalty_min_energy:
                penalty += (penalty_min_energy-real_energy_with_battery[i])*ppa[i]*penalty_multiplier

    fcr = 0.063
    annual_cost = fcr*total_cost + cost_om + penalty

    
    total_energy = np.sum(real_energy_with_battery)

    annual_income = np.sum(real_energy_with_battery*ppa)

    aep = -total_energy
    coe = annual_cost/total_energy
    pro = -(annual_income - annual_cost)/1E6
    
    if objective == "aep":
        obj = aep
    elif objective == "coe":
        obj = coe
    elif objective == "pro":
        obj = pro

    return obj


if __name__=="__main__":

    global rotor_diameter
    global plant
    global wind_wind_spacing
    global wind_solar_spacing
    global sr_array
    global operation
    global boundary_poly
    global wind_capex
    global solar_capex
    global battery_capex
    global om_function
    global HOPP_plant
    global ppa
    global turbine_rating
    global time_array
    global shadow
    global interconnect
    global absolute_min_energy
    global penalty_min_energy
    global penalty_multiplier
    global max_capacity
    global objective


    # Different run considerations
    # 1 interconnect and power requirements: DEFINITELY VARY

    # 3 costs: ...
    # 4 boundary: ...
    # 5 resource: ...
    # 6 PPA: ONLY MATTERS WITH PENALTY AND FOR PROFIT OBJECTIVE

    # 2 turbine design: CONSTANT
    # 7 shadow: CONSTANT


    # interconnect and power requirements
    interconnect = 300.0
    max_capacity = 1E16
    penalty_min_energy = 30.0 # this is the power under which the battery will discharge and a penalty could be applied
    absolute_min_energy = 30.0 # farm must produce at least this amount at all times
    penalty_multiplier = 0.0
    
    # turbine design
    turbine_type = 1
    if turbine_type == 1:
        powercurve_filename = 'high_7r_200d_135h.txt'
        rotor_diameter = 200.0
        hub_height = 135.0
        turbine_rating = 7.0

    HOPP_plant = init_wind_plant(hub_height,rotor_diameter,powercurve_filename)

    # costs
    solar_cost_multiplier = 0.75
    battery_cost_multiplier = 0.5

    capex_cost = np.array([5*1786.0,2*1786.0,1786.0,1622.0,1528.0,1494.0,1470.0,1421.0,1408.0,1400.0]) # $/kW realistic
    capex_size = np.array([0.0,1.0,20.0,50.0,100.0,150.0,200.0,400.0,1000.0,100000.0]) # MW
    cost = capex_size*capex_cost*1000.0
    wind_capex = scipy.interpolate.interp1d(capex_size, cost, kind='cubic')
    solar_capex = scipy.interpolate.interp1d(capex_size, cost*solar_cost_multiplier, kind='cubic')
    battery_capex = scipy.interpolate.interp1d(capex_size, cost*battery_cost_multiplier, kind='cubic')

    def om_function(capacity):
        return 37.0*capacity*1000.0

    # boundary
    bound_type = 1
    if bound_type == 1:
        boundary_poly = Polygon(([0,0],[1000,-1000],[1200,1000],[500,1100],[-200,1000])).buffer(2000)

    outage_start = 3000
    outage_duration = 12
    HOPP_plant = init_wind_plant(hub_height,rotor_diameter,powercurve_filename,outage_start=outage_start,outage_duration=outage_duration)

    np.random.seed(seed=1)
    solar_array = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.15,0.3,0.5,0.7,0.85,0.95,1.00,1.00,0.95,0.85,0.7,0.5,0.3,0.15,0.05,0.0,0.0])
    sr_array = np.array([])
    for i in range(365):
        sr_array = np.append(sr_array,solar_array)

    sr_array = sr_array - np.random.rand(len(sr_array))*0.1
    for i in range(len(sr_array)):
        if sr_array[i] < 0.0:
            sr_array[i] = 0.0

    data = HOPP_plant.site.wind_resource.data["data"]
    nhours = 8760
    time_array = np.arange(nhours)
    for i in range(outage_duration):
        sr_array[outage_start+i] = 0.0


    # PPA
    ppa_type = 1
    if ppa_type == 1:
        ppa = np.zeros(nhours) + 70.0

    elif ppa_type == 2:
        def PPA_curve(hour):
            time = hour%24
            p1 = np.cos((13-time)/(2*np.pi))*70.0
            p2 = 60.0
            return max([p1,p2])
        ppa = np.zeros(nhours)
        for i in range(nhours):
            ppa[i] = PPA_curve(time_array[i])
    
    elif ppa_type == 3:
        def PPA_curve(hour):
            time = hour%24
            p1 = np.cos((time)/(2*np.pi))*70.0
            p2 = np.cos((24-time)/(2*np.pi))*70.0
            p3 = 60.0
            return max([p1,p2,p3])
        ppa = np.zeros(nhours)
        for i in range(nhours):
            ppa[i] = PPA_curve(time_array[i])

    # shadow
    shadow = SimpleShadow()
    shadow.hub_height = hub_height
    shadow.rotor_diameter = rotor_diameter

    shadow_type = 1
    if shadow_type == 1:
        shadow.tower_width = 6.0
        shadow.latitude = 40.966


    operation = "always_stay_charged"

    plant = SimpleHybrid()
    plant.acres_per_MW = 5.0
    D = rotor_diameter
    wind_wind_spacing = D*4.0
    wind_solar_spacing = D*0.5
    
    # coe no min
    
    # wind speed multiplier: 0.6
    # inputs = np.array([9.5238095238095241e+02, 5.0000000000000000e+00,
    #     0.0000000000000000e+00, 7.0000000000000000e+00,
    #     1.2476190476190477e+03, 3.1904761904761904e+03,
    #     2.2439947525641379e-01, 5.8842529067237397e+00,
    #     6.5219769518815974e+02, 2.4862787839403909e+03,
    #     8.7301587301587301e+02, 9.9505141656082242e+02,
    #     2.0991192315633807e+03, 2.0571223282359770e+03,
    #     5.3718387892310229e+03, 4.2857142857142856e+02])
    # opt_val: 71.99168471923504
    # run_time: 4493.82550406456
    # function_calls: 24989

    # wind speed multiplier: 0.8
    # inputs = np.array([ 1.8095238095238096e+03,  1.5000000000000000e+01,
    #         7.0000000000000000e+00,  1.4000000000000000e+01,
    #         3.3523809523809523e+03,  2.4349206349206352e+03,
    #         1.3713301265669733e+00,  4.9866550056980845e+00,
    #         1.6807588593061482e+03, -1.1917369636412104e+03,
    #         6.9841269841269843e+02, -2.9065003858666319e+02,
    #     -5.1420774698144169e+02,  2.0571223282359770e+03,
    #         3.4844359713930958e+03,  3.6190476190476193e+02])
    # opt_val: 65.20042289271726
    # run_time: 8956.88669514656
    # function_calls: 29648

    # wind speed multiplier: 1.0
    # inputs = np.array([ 1.7777777777777778e+03,  2.2000000000000000e+01,
    #         1.3000000000000000e+01,  1.4000000000000000e+01,
    #         8.1587301587301590e+02,  8.1587301587301590e+02,
    #         4.2386567548433729e-01,  2.0943951023931957e+00,
    #     -7.1921719030249164e+02,  1.4215900149036115e+03,
    #         6.9841269841269843e+02,  6.5219769518815974e+02,
    #         2.6011135777257959e+02,  1.2857014551474856e+03,
    #         4.7910994606655067e+03,  3.6190476190476193e+02])
    # opt_val: 58.69435155449355
    # run_time: 9483.828630924225
    # function_calls: 32982

    # wind speed multiplier: 1.2
    # inputs = np.array([ 1.1111111111111111e+03,  2.1000000000000000e+01,
    #         7.0000000000000000e+00,  6.0000000000000000e+00,
    #         9.7777777777777783e+02,  7.0793650793650795e+02,
    #     -6.2333187571226045e-01,  4.9866550056980845e+00,
    #         9.0933798621765709e+02,  8.4085068633809578e+02,
    #         7.4603174603174602e+02,  1.5950454289629824e+03,
    #     -7.0778752316994724e+02,  2.0571223282359770e+03,
    #         2.3229573142620638e+03,  3.6190476190476193e+02])
    # opt_val: 53.0882913748826
    # run_time: 9591.466512203217
    # function_calls: 31583

    # wind speed multiplier: 1.4
    # inputs = np.array([ 0.0000000000000000e+00,  2.2000000000000000e+01,
    #         9.0000000000000000e+00,  1.4000000000000000e+01,
    #         8.1587301587301590e+02,  7.6190476190476193e+02,
    #         3.7399912542735647e-01,  5.6847867064958164e+00,
    #     -8.9064405098882298e+02,  9.3764057443234833e+02,
    #         7.6190476190476193e+02, -3.7636346892982874e+02,
    #         3.2605978886944122e+03,  2.1856924737507256e+03,
    #         2.6133269785448219e+03,  3.6190476190476193e+02])
    # opt_val: 48.498805486310786
    # run_time: 15026.689785718918
    # function_calls: 42353

    # wind speed multiplier: 1.6
    inputs = np.array([ 9.5238095238095241e+01,  2.2000000000000000e+01,
            1.0000000000000000e+01,  1.1000000000000000e+01,
            8.1587301587301590e+02,  7.6190476190476193e+02,
            3.7399912542735647e-01,  1.4959965017094254e+00,
        -1.4049246330478172e+03, -9.0136729935845256e+02,
            5.5555555555555554e+02, -9.7635748133198877e+02,
            1.2280102387151064e+03,  1.0285611641179885e+03,
            2.3229573142620638e+03,  3.6190476190476193e+02])
    # opt_val: 45.89025527443973
    # run_time: 22100.969921588898
    # function_calls: 42984


    objective = "coe"

    plt.figure(figsize=(3,3))
    ax1 = plt.gca()
    turbine_x,turbine_y,solar_polys = get_dvs(inputs)
    plot_hybrid(turbine_x,turbine_y,solar_polys)
    ax1.axis("off")
    ax1.axis("equal")
    plt.subplots_adjust(top=0.99,bottom=0.01,left=0.01,right=0.99)
    plt.savefig("figures/resource1.6.pdf",transparent=True)

    plt.show()



