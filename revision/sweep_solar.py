import time
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from make_solar import PlaceSolar
from plotting_functions import plot_poly, plot_turbines
import matplotlib.pyplot as plt
from gradient_free import GeneticAlgorithm
import scipy.interpolate
from boundary_grid import discrete_grid


def coe_objective(turbine_x, turbine_y, target_solar_capacity, battery_capacity):
    global hybrid_plant
    global wind_function
    global solar_function
    global battery_function
    global turbine_rating
    global wind_om
    global solar_om
    global battery_om

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
    
if __name__=="__main__":
    global hybrid_plant
    global wind_function
    global solar_function
    global battery_function
    global wind_om
    global solar_om
    global battery_om

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

    wind_speed_multiplier = 1.0
    data = site.wind_resource.data["data"]
    nhours = 8760
    for k in range(nhours):
        data[k][2] = data[k][2]*wind_speed_multiplier

    site.wind_resource.data["data"] = data  

    # I'm going to work with the assumption that our plant footprint will be square
    # we can change this of course, but I'll need to alter the solar placement algorithm

    N = 5
    min_turbine_spacing = 5*rotor_diameter
    corner = min_turbine_spacing*(N-1)/2
    
    # change the size if you want, but keep it a square for now so it lines up with the solar grid
    boundary_polygon = Polygon(((-corner,-corner),(-corner,corner),(corner,corner),(corner,-corner)))
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
    min_power = 3000.0
    battery_dispatch = SimpleDispatch(ntime, min_power)

    time_array = np.arange(ntime)

    outage_start = 3008
    outage_duration = 12
    interconnect = 30000.0
    hybrid_plant = HybridPlant(rotor_diameter, hub_height, powercurve_data, place_solar, battery_dispatch, 
                               site, time_array, outage_start, outage_duration, interconnect)


    fcr = 0.063
    wind_cost = np.array([5*1786.0,1786.0,1622.0,1528.0,1494.0,1470.0,1421.0,1408.0])*1555.0/1494.0*fcr # $/kW/year realistic
    wind_capacity = np.array([0.0,20.0,50.0,100.0,150.0,200.0,400.0,1000.0])*1000.0 # MW
    wind_function = scipy.interpolate.interp1d(wind_capacity, wind_cost*wind_capacity, kind='cubic')

    battery_capacity = wind_capacity
    battery_cost = 1284/1555 * wind_cost
    battery_function = scipy.interpolate.interp1d(battery_capacity, battery_cost*battery_capacity, kind='cubic')

    solar_cost_multiplier = 1.2
    solar_capacity = wind_capacity
    solar_cost = 1075/1555 * wind_cost * solar_cost_multiplier
    solar_function = scipy.interpolate.interp1d(solar_capacity, solar_cost*solar_capacity, kind='cubic')

    def wind_om(capacity_kw):
        return 42.0*capacity_kw

    def battery_om(capacity_kw):
        return 32.1*capacity_kw

    def solar_om(capacity_kw):
        return 13.0*capacity_kw

    turbine_x = np.array([ 1037.17647059,   854.56618167,   375.29411765,  1151.22795676,
         671.95589275,   192.68382873,  -286.58823529,   489.34560382,
          10.0735398 ,  -469.19852422,  -948.47058824,  -172.53674912,
        -651.80881314, -1131.08087716])
    turbine_y = np.array([-1153.17647059,  -429.88235294, -1153.17647059,  1016.70588235,
            293.41176471,  -429.88235294, -1153.17647059,  1016.70588235,
            293.41176471,  -429.88235294, -1153.17647059,  1016.70588235,
            293.41176471,  -429.88235294])
    battery_capacity = 228705.8823529412
    # battery_capacity = 0.0
    solar_capacity = 11811.377777777852
    coe = coe_objective(turbine_x, turbine_y, solar_capacity, battery_capacity)
    print(coe)
    # N = 100
    # solar_cap = np.linspace(0,10000,N)
    # coe = np.zeros(N)
    # for i in range(N):
    #     print(i)
    #     coe[i] = coe_objective(turbine_x, turbine_y, solar_cap[i], battery_capacity)

    # import matplotlib.pyplot as plt
    # plt.plot(solar_cap,coe)
    # plt.show()