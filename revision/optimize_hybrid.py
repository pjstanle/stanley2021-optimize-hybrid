import time
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from make_solar import PlaceSolar
from plotting_functions import plot_poly, plot_turbines
import matplotlib.pyplot as plt
from gradient_free import GeneticAlgorithm
import scipy.interpolate
from boundary_grid import discrete_grid


def get_wind_wind_distance(turbine_x, turbine_y):

    #calculate the spacing between each turbine and every other turbine (without repeating)
    nturbs = len(turbine_x)

    if nturbs == 1:
        spacing = np.array([1E6])
    
    else:
        npairs = int((nturbs*(nturbs-1))/2)
        spacing = np.zeros(npairs)

        ind = 0
        for i in range(nturbs):
            for j in range(i,nturbs):
                if i != j:
                    spacing[ind] = np.sqrt((turbine_x[i]-turbine_x[j])**2+(turbine_y[i]-turbine_y[j])**2)
                    ind += 1

    return spacing


def dummy_objective(turbine_x, turbine_y, solar_capacity, battery_storage):
    
    return -(np.sum(turbine_x) + np.sum(turbine_y) + solar_capacity + battery_storage)


def coe_objective(turbine_x, turbine_y, target_solar_capacity, battery_capacity):
    global hybrid_plant
    global wind_function
    global solar_function
    global battery_function
    global turbine_rating
    global wind_om
    global solar_om
    global battery_om

    if target_solar_capacity == 0.0 and len(turbine_x) == 0:
        return 1E20

    hybrid_plant.update_turbine_locs(turbine_x, turbine_y)
    hybrid_plant.solar_capacity = target_solar_capacity+1E-6
    hybrid_plant.battery_size = battery_capacity

    hybrid_plant.evaluate_plant()
    hybrid_power = hybrid_plant.hybrid_power_with_battery
    # print("wind_power: ", np.sum(hybrid_plant.wind_power))
    # print("solar_power: ", np.sum(hybrid_plant.solar_power_with_losses))
    yearly_energy = np.sum(hybrid_power)/1E3

    actual_solar_capacity = hybrid_plant.solar_capacity_kw

    wind_capacity = len(turbine_x) * turbine_rating
    wind_capex = wind_function(wind_capacity)
    solar_capex = solar_function(actual_solar_capacity)
    battery_capex = battery_function(battery_capacity)

    om_wind = wind_om(wind_capacity)
    om_solar = solar_om(actual_solar_capacity)
    om_battery = battery_om(battery_capacity)

    # print("wind_cost: ", wind_capex + om_wind)
    # print("solar_cost: ", solar_capex + om_solar)

    yearly_cost = wind_capex + solar_capex + battery_capex + om_wind + om_solar + om_battery

    coe = yearly_cost / yearly_energy

    return coe

    
def objective_function(inputs):
    global possible_turbine_locs_x
    global possible_turbine_locs_y
    global place_solar
    global solar_cell_area
    global solar_kw_per_km2
    global boundary_polygon
    global min_turbine_spacing

    # nlocs = len(possible_turbine_locs_x)

    # setting up GA to select possible points from predetermined possible turbine locations
    # could also do a grid, boundary-grid, or potentially even fix turbine locations. 
    # BUT, our problem is going to be small, so I think this is the best approach
    # turbine_bool = [bool(y) for y in inputs[0:nlocs]]
    # target_solar_capacity = inputs[nlocs]
    # battery_capacity = inputs[nlocs+1]
    
    # I think these are all the design variables we'll have?

    # turbine_x = possible_turbine_locs_x[turbine_bool]
    # turbine_y = possible_turbine_locs_y[turbine_bool]

    x_spacing = inputs[0]
    y_spacing = inputs[1]
    shear = inputs[2]
    rotation = inputs[3]
    center_x = inputs[4]
    center_y = inputs[5]
    target_solar_capacity = inputs[6]
    battery_capacity = inputs[7]

    turbine_x, turbine_y = discrete_grid(x_spacing,y_spacing,shear,rotation,center_x,center_y,0.0,boundary_polygon)

    if len(turbine_x) > 1:
        spacing = get_wind_wind_distance(turbine_x, turbine_y)

        if min(spacing) < min_turbine_spacing:
            return 1E16

    # GA is setup to minimize
    # not sure if you'll have it setup with solar capacity or the actual solar polygons (like if we do shadow/flicker)
    # objective = dummy_objective(turbine_x, turbine_y, solar_capacity_kw, battery_capacity)
    objective = coe_objective(turbine_x, turbine_y, target_solar_capacity, battery_capacity)

    # plt.cla()
    # plot_poly(solar_polygons,ax=plt.gca())
    # plot_turbines(turbine_x,turbine_y,rotor_diameter/2,ax=plt.gca())
    # plt.pause(0.001)

    return objective




if __name__=="__main__":

    global possible_turbine_locs_x
    global possible_turbine_locs_y
    global place_solar
    global solar_cell_area
    global solar_kw_per_km2
    global wind_function
    global solar_function
    global battery_function
    global wind_om
    global solar_om
    global battery_om
    global turbine_rating
    global min_turbine_spacing

    outage_start = 3000
    outage_duration = 12

    # setting up the power calculations
    import os
    from dotenv import load_dotenv
    from hybrid.sites import SiteInfo
    from hybrid.sites import flatirons_site as sample_site
    from hybrid.keys import set_developer_nrel_gov_key
    import json
    import matplotlib.pyplot as plt
    from init_plant import init_hybrid_plant
    from simple_dispatch import SimpleDispatch
    from hybrid_plant import HybridPlant

    # Set API key
    load_dotenv()
    NREL_API_KEY = os.getenv("NREL_API_KEY")
    set_developer_nrel_gov_key(NREL_API_KEY)  # Set this key manually here if you are not setting it using the .env
    
    # big turbine
    # powercurve_filename = 'turbine_data/high_7r_200d_135h.txt'
    # rotor_diameter = 200.0
    # hub_height = 135.0
    # smaller turbine
    powercurve_filename = 'turbine_data/low_2_43r_116d_88h.txt'
    rotor_diameter = 116.0
    hub_height = 88.0
    turbine_rating = 2430.0

    # Get resource
    lat = 40.966
    lon = -84.598
    year = 2012

    sample_site['year'] = year
    sample_site['lat'] = lat
    sample_site['lon'] = lon

    # Import powercurve
    powercurve_file = open(powercurve_filename)
    powercurve_data = json.load(powercurve_file)
    powercurve_file.close()

    site = SiteInfo(sample_site, hub_height=hub_height)


    # I'm going to work with the assumption that our plant footprint will be square
    # we can change this of course, but I'll need to alter the solar placement algorithm

    N = 5
    min_turbine_spacing = 5*rotor_diameter
    corner = min_turbine_spacing*(N-1)/2
    
    # change the size if you want, but keep it a square for now so it lines up with the solar grid
    boundary_polygon = Polygon(((-corner,-corner),(-corner,corner),(corner,corner),(corner,-corner)))
    # tx = np.arange(N)*min_turbine_spacing - corner
    # ty = np.arange(N)*min_turbine_spacing - corner
    # possible_turbine_locs_x, possible_turbine_locs_y = np.meshgrid(tx, ty)
    # possible_turbine_locs_x = np.ndarray.flatten(possible_turbine_locs_x)
    # possible_turbine_locs_y = np.ndarray.flatten(possible_turbine_locs_y)

    # if you want to see
    # plot_poly(boundary_polygon,ax=plt.gca())
    # plot_turbines(possible_turbine_locs_x,possible_turbine_locs_y,rotor_diameter/2,ax=plt.gca())
    # plt.axis("equal")
    # plt.show()

    # setting up the solar placement algorithm
    ngrid = 25 # dunno what this should be, affects how coarse the solar grid is
    x = np.linspace(-corner,corner,ngrid)
    y = np.linspace(-corner,corner,ngrid)
    xlocs = [[i for i in x] for j in y]
    ylocs = [[j for i in x] for j in x]
    grid_locs = np.zeros((np.shape(xlocs)[0],np.shape(xlocs)[1],2))
    grid_locs[:,:,0] = xlocs[:]
    grid_locs[:,:,1] = ylocs[:]
    grid_locs = np.ndarray.tolist(grid_locs)

    min_solar_spacing = 2*rotor_diameter
    place_solar = PlaceSolar(grid_locs,min_solar_spacing)

    dx = x[1] - x[0]
    dy = y[1] - y[0]
    solar_cell_area = dx * dy

    solar_kw_per_km2 = 4000.0
    
    ntime = 8760
    # battery
    min_power = 0.0
    battery_dispatch = SimpleDispatch(ntime, min_power)

    time_array = np.arange(ntime)
    hybrid_plant = HybridPlant(rotor_diameter, hub_height, powercurve_data, place_solar, battery_dispatch, site, time_array, outage_start, outage_duration)


    fcr = 0.063
    wind_cost = np.array([5*1786.0,1786.0,1622.0,1528.0,1494.0,1470.0,1421.0,1408.0])*1555.0/1494.0*fcr # $/kW/year realistic
    wind_capacity = np.array([0.0,20.0,50.0,100.0,150.0,200.0,400.0,1000.0])*1000.0 # MW
    wind_function = scipy.interpolate.interp1d(wind_capacity, wind_cost*wind_capacity, kind='cubic')

    battery_capacity = wind_capacity
    battery_cost = 1284/1555 * wind_cost
    battery_function = scipy.interpolate.interp1d(battery_capacity, battery_cost*battery_capacity, kind='cubic')

    solar_cost_multiplier = 0.5
    solar_capacity = wind_capacity
    solar_cost = 1075/1555 * wind_cost * solar_cost_multiplier
    solar_function = scipy.interpolate.interp1d(solar_capacity, solar_cost*solar_capacity, kind='cubic')

    # x = np.linspace(0,1000,100)*1000
    # plt.plot(x,wind_function(x)/x,label="wind")
    # plt.plot(x,battery_function(x)/x,label="battery")
    # plt.plot(x,solar_function(x)/x,label="solar")
    # plt.legend()
    # plt.show()

    def wind_om(capacity_kw):
        return 42.0*capacity_kw

    def battery_om(capacity_kw):
        return 32.1*capacity_kw

    def solar_om(capacity_kw):
        return 13.0*capacity_kw

    # test objective function
    # turbs = np.random.randint(0,high=2,size=len(possible_turbine_locs_x))
    # turbs = np.zeros(len(possible_turbine_locs_x),dtype=int)
    # turbs[0:3] = 1

    # x_spacing = min_turbine_spacing*10
    # y_spacing = min_turbine_spacing*10
    # shear = 0.0
    # rotation = np.pi/6
    # center_x = 1000.0
    # center_y = 2000.0
    # target_solar_capacity = 0.0
    # battery_capacity = 0.0

    # x = np.array([x_spacing, y_spacing, shear, rotation, center_x, center_y, 
    #               target_solar_capacity, battery_capacity])

    # # x = np.append(turbs,solar_capacity)
    # # x = np.append(x,battery_capacity)
    # coe = objective_function(x)
    # print(coe)
    # # turbine_bool = [bool(y) for y in turbs]
    # # turbine_x = possible_turbine_locs_x[turbine_bool]   
    # # turbine_y = possible_turbine_locs_y[turbine_bool]
    # turbine_x, turbine_y = discrete_grid(x_spacing,y_spacing,shear,rotation,center_x,center_y,0.0,boundary_polygon)
    # turbine_locs = np.zeros((len(turbine_x),2))
    # turbine_locs[:,0] = turbine_x
    # turbine_locs[:,1] = turbine_y
    # hybrid_plant.solar.set_turbine_locs(turbine_locs)
    # nsolar_cells = int(np.floor(target_solar_capacity/(solar_kw_per_km2*solar_cell_area/1E6)))
    # hybrid_plant.solar.nsolar_cells = nsolar_cells
    # hybrid_plant.solar.place_solar()
    # solar_polygons = hybrid_plant.solar.solar_geometry
    # plot_poly(boundary_polygon,ax=plt.gca(),fill=False)
    # plot_poly(solar_polygons,ax=plt.gca())
    # plot_turbines(turbine_x,turbine_y,rotor_diameter/2,ax=plt.gca())
    # plt.axis("equal")
    # plt.show()






    # # setting up the optimizer

    # nlocs = len(possible_turbine_locs_x)
    # # bits per design variable. Turbine locs are Bool so they just need 1 each. The other 2 are floats so...8 bits each should be good
    # bits = np.ones(nlocs+2, dtype=int)
    # bits[-1] = 8
    # bits[-2] = 8

    # # bounds for each variable. 0 to 2 (exclusive) for the turbine locs. 
    # bounds = np.zeros((nlocs+2, 2), dtype=int)
    # bounds[:, 1] = 2
    # bounds[-1,0] = 0.0
    # bounds[-1,1] = 25000.0
    # bounds[-2,0] = 0.0 
    # bounds[-2,1] = 25000.0

    # # variable type for each variable (int or float)
    # variable_type = np.array([])
    # for _ in range(nlocs):
    #     variable_type = np.append(variable_type, "int")

    # variable_type = np.append(variable_type, "float")
    # variable_type = np.append(variable_type, "float")

    # x_spacing = min_turbine_spacing*10
    # y_spacing = min_turbine_spacing*10
    # shear = 0.0
    # rotation = np.pi/6
    # center_x = 1000.0
    # center_y = 2000.0
    # target_solar_capacity = 0.0
    # battery_capacity = 0.0
    area = (2*corner)**2
    max_solar = area/1E6*solar_kw_per_km2
    bits = np.array([8,8,8,8,8,8,8,8])
    bounds = np.array([(min_turbine_spacing,10*min_turbine_spacing),(min_turbine_spacing,10*min_turbine_spacing),
                        (0,np.pi/4),(-np.pi,np.pi),(-1.5*corner,1.5*corner),(-1.5*corner,1.5*corner),
                        (0,max_solar*1.1),(0,10*turbine_rating*25)])
    variable_type = np.array(["float","float","float","float","float","float","float","float"])

    ga = GeneticAlgorithm()
    ga.bits = bits
    ga.bounds = bounds
    ga.variable_type = variable_type
    ga.objective_function = objective_function

    # probably worth playing with these parameters
    ga.max_generation = 500
    ga.population_size = 10
    ga.crossover_rate = 0.1
    ga.mutation_rate = 0.02
    ga.tol = 0.0
    ga.convergence_iters = 20

    ga.optimize_ga(print_progress=True)

    # unpack results
    solution_history = ga.solution_history
    function_value = ga.optimized_function_value
    optimal_design_variables = ga.optimized_design_variables

    # get optimal design variables into meaningful values
    # turbine_bool = [bool(y) for y in optimal_design_variables[0:nlocs]]
    # target_solar_capacity = optimal_design_variables[nlocs]
    # battery_capacity = optimal_design_variables[nlocs+1]
    x_spacing = optimal_design_variables[0]
    y_spacing = optimal_design_variables[1]
    shear = optimal_design_variables[2]
    rotation = optimal_design_variables[3]
    center_x = optimal_design_variables[4]
    center_y = optimal_design_variables[5]
    target_solar_capacity = optimal_design_variables[6]
    battery_capacity = optimal_design_variables[7]

    # turbine_bool = [bool(y) for y in optimal_design_variables[0:nlocs]]
    # target_solar_capacity = optimal_design_variables[nlocs]
    # battery_capacity = optimal_design_variables[nlocs+1]
    # turbine_x = possible_turbine_locs_x[turbine_bool]
    # turbine_y = possible_turbine_locs_y[turbine_bool]

    turbine_x, turbine_y = discrete_grid(x_spacing,y_spacing,shear,rotation,center_x,center_y,0.0,boundary_polygon)

    turbine_locs = np.zeros((len(turbine_x),2))
    turbine_locs[:,0] = turbine_x
    turbine_locs[:,1] = turbine_y
    place_solar.set_turbine_locs(turbine_locs)

    nsolar_cells = int(np.floor(target_solar_capacity/(solar_kw_per_km2*solar_cell_area/1E6)))
    place_solar.nsolar_cells = nsolar_cells
    place_solar.place_solar()
    solar_polygons = place_solar.solar_geometry
    solar_capacity_kw = solar_polygons.area/1E6*solar_kw_per_km2

    # print("turbine_bool: ", turbine_bool)
    print("turbine_x: ", repr(turbine_x))
    print("turbine_y: ", repr(turbine_y))
    print("wind_capacity: ", len(turbine_x)*turbine_rating)
    print("solar_capacity: ", solar_capacity_kw)
    print("battery_capacity: ", battery_capacity)

    plt.figure(1)
    plt.plot(solution_history)
    plt.xlabel("generation")
    plt.ylabel("function value")

    plt.figure(2)
    plot_poly(solar_polygons,ax=plt.gca())
    plot_turbines(turbine_x,turbine_y,rotor_diameter/2,ax=plt.gca())
    plt.axis("equal")

    plt.show()