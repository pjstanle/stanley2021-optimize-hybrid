
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
    # plt.plot(real_energy_with_battery)
    # plt.show()
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
    ppa_type = 3
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

    # 0 hours
    inputs0 = np.array([ 5.3968253968253975e+02,  2.2000000000000000e+01,  1.3000000000000000e+01,
        1.1000000000000000e+01 , 9.7777777777777783e+02  ,9.7777777777777783e+02,
        4.7373222554131811e-01 , 3.8895909044445061e+00  ,6.5219769518815974e+02,
        3.0670181125059071e+03 , 8.2539682539682542e+02  ,1.7664722896493140e+03,
        -1.2704819460443105e+02,  6.4285072757374282e+03 , 1.4518483214137900e+03,
        3.0476190476190476e+02])
    
    # 3 hours
    inputs3 = np.array([ 6.0317460317460313e+02,  2.2000000000000000e+01,
        1.3000000000000000e+01,  9.0000000000000000e+00,
        8.1587301587301590e+02,  1.1396825396825398e+03,
        1.7453292519943298e-01,  5.1861212059260078e+00,
       -3.3509747557166065e+01,  5.5048102205533769e+02,
        8.0952380952380952e+02, -3.7636346892982874e+02,
        1.6151697910921166e+03,  1.5428417461769827e+03,
        5.8073932856551601e+03,  3.0476190476190476e+02])
    
    # 6 hours
    inputs6 = np.array([ 3.1746031746031747e+01,  2.1000000000000000e+01 , 1.1000000000000000e+01,
        1.1000000000000000e+01 , 1.5174603174603176e+03 , 6.0000000000000000e+02,
        -8.2279807594018384e-01,  4.4879895051282759e+00, -9.7635748133198877e+02,
        1.4215900149036115e+03 , 7.4603174603174602e+02 , 1.3379051379334855e+03,
        8.4085068633809578e+02 , 2.3142626192654743e+03 , 2.9036966428275800e+03,
        3.0476190476190476e+02])

    # 9 hours
    inputs9 = np.array([4.7619047619047620e+02, 2.2000000000000000e+01,
       1.0000000000000000e+01, 7.0000000000000000e+00,
       9.2380952380952385e+02, 1.2476190476190477e+03,
       6.2333187571226079e-01, 4.7871888054701612e+00,
       9.9505141656082242e+02, 2.6011135777257959e+02,
       7.7777777777777783e+02, 3.0934397381549707e+02,
       3.1638080006001601e+03, 1.7999820372064798e+03,
       7.5496112713517077e+03, 3.0476190476190476e+02])
    
    # 12 hours
    inputs12 = np.array([ 6.3492063492063494e+02,  2.1000000000000000e+01  ,5.0000000000000000e+00,
        1.3000000000000000e+01 , 1.1396825396825398e+03,  6.5396825396825398e+02,
        -6.2333187571226045e-01 , 2.5930606029630039e+00, -2.0049186454499772e+03,
        1.5183799029978645e+03 , 5.3968253968253975e+02,  8.2362455587449131e+02,
        1.2280102387151064e+03 , 2.6999730558097199e+03 , 2.3229573142620638e+03,
        3.6190476190476193e+02])

    # 15 hours
    inputs15 = np.array([4.1269841269841271e+02, 2.2000000000000000e+01,
       1.3000000000000000e+01, 1.0000000000000000e+01,
       8.6984126984126988e+02, 7.6190476190476193e+02,
       6.2333187571226079e-01, 3.0917261035328125e+00,
       2.8807468841104678e+03, 2.5830686720346439e+03,
       7.6190476190476193e+02, 1.2521917075903198e+03,
       3.5690124586683214e+02, 2.6999730558097199e+03,
       2.3229573142620638e+03, 4.5714285714285711e+02])

    # # 18 hours

    # inputs18 = np.array([ 8.8888888888888891e+02,  2.2000000000000000e+01,
    #     1.5000000000000000e+01,  1.1000000000000000e+01,
    #     1.4634920634920636e+03,  7.0793650793650795e+02,
    #     5.2359877559829915e-01,  3.6901247042165828e+00,
    #     6.5219769518815974e+02, -3.0258306510178500e+01,
    #     7.6190476190476193e+02,  6.5219769518815974e+02,
    #    -7.0778752316994724e+02,  1.6714118916917314e+03,
    #     6.2429477820792963e+03,  5.4285714285714289e+02])

    # # 21 hours
    # inputs21 = np.array([ 6.3492063492063494e+02,  5.0000000000000000e+01,
    #     1.2000000000000000e+01,  5.0000000000000000e+00,
    #     3.0825396825396824e+03,  2.2190476190476193e+03,
    #     1.5707963267948966e+00,  1.9946620022792338e-01,
    #     2.9664603144536336e+03, -3.0258306510178500e+01,
    #     2.3809523809523810e+02,  5.6648426484499396e+02,
    #    -2.2383808269868359e+02,  5.9142266936784335e+03,
    #     6.9688719427861915e+03,  3.9365079365079362e+02])


    # # 24 hours
    # inputs24 = np.array([ 1.4603174603174602e+03,  2.2000000000000000e+01,
    #     1.5000000000000000e+01,  1.2000000000000000e+01,
    #     8.1587301587301590e+02,  7.6190476190476193e+02,
    #    -4.7373222554131789e-01,  4.9866550056980845e+00,
    #    -2.9065003858666319e+02, -2.0628459564894847e+03,
    #     7.7777777777777783e+02,  1.5950454289629824e+03,
    #     6.4727091014959024e+02,  2.0571223282359770e+03,
    #     3.1940663071103377e+03,  9.4285714285714289e+02])


    # inputs04 = np.array([ 0.0000000000000000e+00,  2.2000000000000000e+01,
    #     7.0000000000000000e+00,  8.0000000000000000e+00,
    #     8.6984126984126988e+02,  1.4095238095238096e+03,
    #     1.0721308262250884e+00,  2.0943951023931957e+00,
    #    -1.1922317790033162e+02,  1.3248001268093594e+03,
    #     8.7301587301587301e+02,  5.6648426484499396e+02,
    #     5.5048102205533769e+02,  4.6285252385309486e+03,
    #     1.4518483214137900e+03,  3.6190476190476193e+02])
    # inputs06 = np.array([ 7.6190476190476193e+02,  2.2000000000000000e+01,
    #     1.5000000000000000e+01,  1.1000000000000000e+01,
    #     1.0857142857142858e+03,  9.2380952380952385e+02,
    #    -3.2413257537037543e-01,  3.9893240045584677e+00,
    #    -5.4779032961616031e+02, -5.1420774698144169e+02,
    #     6.6666666666666674e+02,  3.1378871751399652e+03,
    #    -2.0628459564894847e+03,  4.4999550930161995e+03,
    #     4.9362842928068858e+03,  3.6190476190476193e+02])
    # inputs08 = np.array([ 2.5396825396825398e+02,  2.2000000000000000e+01,
    #     9.0000000000000000e+00,  8.0000000000000000e+00,
    #     9.2380952380952385e+02,  8.6984126984126988e+02,
    #     8.7266462599716510e-01,  5.7845198066097785e+00,
    #    -1.2334977723614859e+03, -5.1420774698144169e+02,
    #     7.3015873015873012e+02,  9.9505141656082242e+02,
    #    -4.1741785888718914e+02,  1.4142716006622341e+03,
    #     4.2103601320999906e+03,  3.6190476190476193e+02])
    # inputs10 = np.array([ 6.6666666666666674e+02,  2.2000000000000000e+01,
    #     1.5000000000000000e+01,  1.5000000000000000e+01,
    #     1.1396825396825398e+03,  8.1587301587301590e+02,
    #     1.1219973762820690e+00,  3.6901247042165828e+00,
    #    -1.8334917847636459e+03, -1.9660560683952319e+03,
    #     6.1904761904761904e+02,  1.5950454289629824e+03,
    #     6.4727091014959024e+02,  2.0571223282359770e+03,
    #     2.9036966428275800e+03,  3.6190476190476193e+02])
    # inputs12 = np.array([ 1.5873015873015875e+03,  2.2000000000000000e+01,
    #     1.5000000000000000e+01,  1.4000000000000000e+01,
    #     8.1587301587301590e+02,  8.1587301587301590e+02,
    #    -1.2466637514245194e-01,  4.6874557053561992e+00,
    #     1.4236185682766513e+03,  1.7119596791863696e+03,
    #     7.6190476190476193e+02, -1.1922317790033162e+02,
    #    -4.1741785888718914e+02,  2.3142626192654743e+03,
    #     2.1777724821206848e+03,  3.6190476190476193e+02])
    # inputs14 = np.array([1.9047619047619048e+02, 2.2000000000000000e+01,
    #    1.5000000000000000e+01, 9.0000000000000000e+00,
    #    1.0317460317460318e+03, 6.5396825396825398e+02,
    #    6.7319842576924138e-01, 3.8895909044445061e+00,
    #    6.5219769518815974e+02, 1.5183799029978645e+03,
    #    7.7777777777777783e+02, 1.3379051379334855e+03,
    #    1.3248001268093594e+03, 2.0571223282359770e+03,
    #    2.3229573142620638e+03, 3.6190476190476193e+02])
    # inputs16 = np.array([ 1.6507936507936508e+03,  2.2000000000000000e+01,
    #     1.4000000000000000e+01,  1.5000000000000000e+01,
    #     8.1587301587301590e+02,  8.1587301587301590e+02,
    #     2.4933275028490520e-02,  6.2831853071795862e+00,
    #     5.6648426484499396e+02,  6.6531581584074047e+01,
    #     7.9365079365079373e+02, -1.1922317790033162e+02,
    #    -4.1741785888718914e+02,  2.1856924737507256e+03,
    #     1.1614786571310319e+03,  3.6190476190476193e+02])


    # inputs0 = np.array([ 8.5714285714285711e+02,  0.0000000000000000e+00,
    #         8.0000000000000000e+00,  0.0000000000000000e+00,
    #         2.5968253968253966e+03,  1.0857142857142858e+03,
    #         1.5209297767379160e+00,  6.2831853071795862e+00,
    #         5.6648426484499396e+02,  2.5830686720346439e+03,
    #         8.4126984126984132e+02, -2.0493660824349740e+02,
    #     -2.2383808269868359e+02,  1.9285521827212283e+03,
    #         3.1940663071103377e+03,  0.0000000000000000e+00])
    # inputs10 = np.array([ 1.6825396825396826e+03,  1.5000000000000000e+01,
    #         7.0000000000000000e+00,  2.0000000000000000e+00,
    #         2.1650793650793648e+03,  1.4095238095238096e+03,
    #     -2.7426602531339461e-01,  2.7925268031909276e+00,
    #     -3.3509747557166065e+01,  9.3764057443234833e+02,
    #         1.0000000000000000e+03, -2.9065003858666319e+02,
    #         7.4406079824384278e+02,  1.6714118916917314e+03,
    #         3.6296208035344748e+03,  1.2698412698412699e+02])
    # inputs20 = np.array([ 3.1746031746031747e+01,  2.2000000000000000e+01,
    #         7.0000000000000000e+00,  1.1000000000000000e+01,
    #         1.5714285714285716e+03,  1.0317460317460318e+03,
    #     -2.2439947525641379e-01,  5.9839860068377015e-01,
    #     -1.3192112027046517e+03, -2.1596358445837373e+03,
    #         7.6190476190476193e+02,  2.2363054347233128e+02,
    #     -4.1741785888718914e+02,  2.9571133468392168e+03,
    #         2.0325876499793058e+03,  2.5396825396825398e+02])
    # inputs30 = np.array([ 6.3492063492063494e+02 , 2.1000000000000000e+01,  5.0000000000000000e+00,
    #     1.3000000000000000e+01 , 1.1396825396825398e+03,  6.5396825396825398e+02,
    #     -6.2333187571226045e-01,  2.5930606029630039e+00, -2.0049186454499772e+03,
    #     1.5183799029978645e+03 , 5.3968253968253975e+02 , 8.2362455587449131e+02,
    #     1.2280102387151064e+03 , 2.6999730558097199e+03 , 2.3229573142620638e+03,
    #     3.6190476190476193e+02])
    # inputs40 = np.array([ 1.0158730158730159e+03,  2.2000000000000000e+01,
    #         9.0000000000000000e+00,  1.1000000000000000e+01,
    #         1.0317460317460318e+03,  8.1587301587301590e+02,
    #         7.2306497582622242e-01,  6.2831853071795862e+00,
    #     -4.6207689927299452e+02, -1.2885268517354632e+03,
    #         7.3015873015873012e+02,  7.3791112553132552e+02,
    #         6.4727091014959024e+02,  4.1142446564719539e+03,
    #         1.5970331535551688e+03,  5.0793650793650795e+02])
    # inputs50 = np.array([ 9.8412698412698410e+02,  2.2000000000000000e+01,
    #         1.0000000000000000e+01,  1.3000000000000000e+01,
    #         1.0857142857142858e+03,  6.0000000000000000e+02,
    #     -7.2306497582622220e-01,  0.0000000000000000e+00,
    #     -4.6207689927299452e+02,  6.6531581584074047e+01,
    #         7.7777777777777783e+02,  9.9505141656082242e+02,
    #         6.6531581584074047e+01,  2.8285432013244681e+03,
    #         2.3229573142620638e+03,  6.0317460317460313e+02])

    objective = "coe"

    inputs = np.array([inputs0,inputs3,inputs6,inputs9,inputs12,inputs15])
    # inputs = np.array([inputs0,inputs3,inputs6,inputs9,inputs12,inputs15,inputs18,inputs21,inputs24])
    # inputs = np.array([inputs04,inputs06,inputs08,inputs10,inputs12,inputs14,inputs16])
    # inputs = np.array([inputs0,inputs10,inputs20,inputs30,inputs40,inputs50])
    n = 6
    solar_capacity = np.zeros(n)
    wind_capacity = np.zeros(n)
    battery_capacity = np.zeros(n)
    coe = np.zeros(n)
    outage = np.array([0,3,6,9,12,15])
    # price = np.array([0.4,0.6,0.8,1.0,1.2,1.4,1.6])
    # min_power = np.array([0,10,20,30,40,50])
    smpa = 4046.86 # square meters per acre

    # plt.figure(1,figsize=(12,2))
    for i in range(n):
        # solar_cost_multiplier = price[i]
        # solar_capex = scipy.interpolate.interp1d(capex_size, cost*solar_cost_multiplier, kind='cubic')

        # absolute_min_energy = min_power[i]
        # penalty_min_energy = min_power[i]
        outage_duration = outage[i]
        # outage_duration = 12
        HOPP_plant = init_wind_plant(hub_height,rotor_diameter,powercurve_filename,outage_start=outage_start,outage_duration=outage_duration)

        np.random.seed(seed=1)
        solar_array = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.15,0.3,0.5,0.7,0.85,0.95,1.00,1.00,0.95,0.85,0.7,0.5,0.3,0.15,0.05,0.0,0.0])
        sr_array = np.array([])
        for k in range(365):
            sr_array = np.append(sr_array,solar_array)

        sr_array = sr_array - np.random.rand(len(sr_array))*0.1
        for k in range(len(sr_array)):
            if sr_array[k] < 0.0:
                sr_array[k] = 0.0

        turbine_x,turbine_y,solar_polys = get_dvs(inputs[i,:])
        solar_capacity[i] = solar_polys.area/smpa/plant.acres_per_MW
        wind_capacity[i] = len(turbine_x)*turbine_rating
        battery_capacity[i] = inputs[i][-1]
        coe[i] = objective_function(inputs[i])
        # plt.subplot(1,n,i+1)
        # # plt.title("%s"%price[i])
        # plot_hybrid(turbine_x,turbine_y,solar_polys)
        # plt.axis("equal")
        # plt.axis("off")

    # plt.subplots_adjust(left=0,right=1.0,top=1.0,bottom=0)
    # plt.show()
    
    print("solar capacity: ", repr(solar_capacity))
    print("wind capacity: ", repr(wind_capacity))
    print("battery capacity: ", repr(battery_capacity))
    print("coe: ", repr(coe))

    xaxis = outage
    plt.figure(1)
    plt.plot(xaxis,solar_capacity)
    plt.title("solar capacity")

    plt.figure(2)
    plt.plot(xaxis,wind_capacity)
    plt.title("wind capacity")

    plt.figure(3)
    plt.plot(xaxis,battery_capacity)
    plt.title("battery capacity")

    plt.figure(4)
    plt.plot(xaxis,coe)
    plt.title("coe")

    plt.figure(5)
    plt.plot(xaxis,solar_capacity+wind_capacity)
    plt.title("total capacity")

    plt.show()

    
