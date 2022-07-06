import time
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from make_solar import PlaceSolar
from plotting_functions import plot_poly, plot_turbines
import matplotlib.pyplot as plt
from gradient_free import GeneticAlgorithm
import scipy.interpolate
from boundary_grid import discrete_grid
import matplotlib.pyplot as plt


def coe_objective(turbine_x, turbine_y, target_solar_capacity, battery_capacity):
    global hybrid_plant
    global wind_function
    global solar_function
    global battery_function
    global turbine_rating
    global wind_om
    global solar_om
    global battery_om
    global min_power

    if target_solar_capacity == 0.0 and len(turbine_x) == 0:
        return 1E20

    hybrid_plant.update_turbine_locs(turbine_x, turbine_y)
    hybrid_plant.solar_capacity = target_solar_capacity+1E-6
    hybrid_plant.battery_size = battery_capacity

    hybrid_plant.evaluate_plant()
    hybrid_power = hybrid_plant.hybrid_power_with_battery

    tol = 1E-6
    if np.min(hybrid_power) < min_power-tol:
        # return 1E20
        print("FAIL")
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

    yearly_cost = wind_capex + solar_capex + battery_capex + om_wind + om_solar + om_battery

    coe = yearly_cost / yearly_energy

    return coe, yearly_cost, yearly_energy
    

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

    outage_start = 3008
    outage_duration = 12
    interconnect = 100000.0
    min_power = 12000.0
    wind_speed_multiplier = 1.0

    powercurve_filename = 'turbine_data/high_7r_200d_135h.txt'
    rotor_diameter = 200.0
    hub_height = 135.0
    turbine_rating = 7000.0

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

    data = site.wind_resource.data["data"]
    nhours = 8760
    for k in range(nhours):
        data[k][2] = data[k][2]*wind_speed_multiplier

    site.wind_resource.data["data"] = data  
    
    # I'm going to work with the assumption that our plant footprint will be square
    # we can change this of course, but I'll need to alter the solar placement algorithm

    N = 5
    min_turbine_spacing = 5*rotor_diameter
    # corner = min_turbine_spacing*(N-1)/2
    corner = 2000.0
    
    # change the size if you want, but keep it a square for now so it lines up with the solar grid
    boundary_polygon = Polygon(((-corner,-corner),(-corner,corner),(corner,corner),(corner,-corner)))

    # setting up the solar placement algorithm
    ngrid = 20 # dunno what this should be, affects how coarse the solar grid is
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

    solar_kw_per_km2 = 1000.0/5.0 * 247.105
    
    ntime = 8760
    # battery
    
    battery_dispatch = SimpleDispatch(ntime, min_power)

    time_array = np.arange(ntime)
    hybrid_plant = HybridPlant(rotor_diameter, hub_height, powercurve_data, place_solar, battery_dispatch, site, 
                            time_array, outage_start, outage_duration, interconnect)

    fcr = 0.063
    wind_cost = np.array([5*1786.0,1786.0,1622.0,1528.0,1494.0,1470.0,1421.0,1408.0])*1555.0/1494.0*fcr # $/kW/year realistic
    wind_capacity = np.array([0.0,20.0,50.0,100.0,150.0,200.0,400.0,1000.0])*1000.0 # MW
    wind_function = scipy.interpolate.interp1d(wind_capacity, wind_cost*wind_capacity, kind='cubic')

    battery_capacity = wind_capacity
    battery_cost = 1284/1555 * wind_cost
    battery_function = scipy.interpolate.interp1d(battery_capacity, battery_cost*battery_capacity, kind='cubic')

    solar_cost_multiplier = 1.0
    solar_capacity = wind_capacity
    solar_cost = 1075/1555 * wind_cost * solar_cost_multiplier
    solar_function = scipy.interpolate.interp1d(solar_capacity, solar_cost*solar_capacity, kind='cubic')

    def wind_om(capacity_kw):
        return 42.0*capacity_kw

    def battery_om(capacity_kw):
        return 32.1*capacity_kw

    def solar_om(capacity_kw):
        return 13.0*capacity_kw

    
    solar1 = 359226.3268698072
    battery1 = 250244.37927663734
    x1 = [-1953.98043925,  -667.4664101 ,   619.04761905,  1905.56164819,
       -1925.48822012,  -638.97419097,   647.53983818,  1934.05386732,
       -1896.99600099,  -610.48197184,   676.03205731,  1962.54608645,
       -1868.50378186,  -581.98975271,   704.52427644,  1991.03830558]
    y1 = [ 1508.26914618,  1492.22981119,  1476.19047619,  1460.15114119,
         365.7672239 ,   349.7278889 ,   333.6885539 ,   317.64921891,
        -776.73469839,  -792.77403339,  -808.81336838,  -824.85270338,
       -1919.23662068, -1935.27595567, -1951.31529067, -1967.35462567]

    solar = np.array([297895.0027700839,367987.9445983399,365797.5401662054,481888.9750692545,300085.40720221744])
    battery = np.array([125122.18963831867,172043.01075268816,250244.37927663734,375366.568914956,688172.0430107526])


    turbine_x_array = np.array([[-1476.19047619, -1645.89965972, -1815.60884325, -1985.31802678,
        -281.59664052,  -451.30582405,  -621.01500758,  -790.72419111,
         743.28801163,   573.5788281 ,   403.86964457,  1937.8818473 ,
        1768.17266377,  1598.46348024,  1428.75429671],[-1813.81121154, -1800.97993974, -1203.30746228, -1190.47619048,
        -592.80371302,     4.86876444,  -579.97244121,    17.70003625,
         615.37251371,  1213.04499117,   628.20378551,  1225.87626297,
        1823.54874043,  1836.38001224],[-1453.44340311,  -411.60773749,   630.22792812,  1672.06359374,
       -1043.96056707,    -2.12490145,  1039.71076416, -1676.31339665,
        -634.47773103,   407.35793458,  1449.1936002 , -1266.83056061,
        -224.99489499,   816.84077062,  1858.67643624],[ 1856.84837624,   861.87633561,  -133.09570501, -1128.06774563,
        1899.58670493,   904.61466431,   -90.35737631, -1085.32941693,
        1942.32503362,   947.352993  ,   -47.61904762],[ 1951.17725614,  1726.43757991,  1278.70885354,   830.98012717,
        1053.96917731,   606.24045094,   158.51172457,  -289.2170018 ,
         -66.22795166,  -513.95667803,  -961.6854044 , -1409.41413077,
       -1857.14285714, -1186.42508063, -1634.153807  ]])

    turbine_y_array = np.array([[-1571.42857143,  -445.48170767,   680.4651561 ,  1806.41201986,
       -1921.05547906,  -795.1086153 ,   330.83824847,  1456.78511223,
       -1144.73552293,   -18.78865916,  1107.1582046 , -1494.36243056,
        -368.41556679,   757.53129697,  1883.47816073],[ 1870.41982854,  -253.2655154 ,  1028.44724871, -1095.23809524,
         186.47466887,  1468.18743299, -1937.21067507,  -655.49791096,
         626.21485315,  1907.92761726, -1497.4704908 ,  -215.75772669,
        1065.95503743, -1057.73030652],[ 1820.21900943,  1859.20174965,  1898.18448988,  1937.1672301 ,
         601.45527865,   640.43801887,   679.42075909,  -656.29119236,
        -617.30845214,  -578.32571192,  -539.3429717 , -1875.05492315,
       -1836.07218292, -1797.0894427 , -1758.10670248],[-1903.28899487, -1689.20690467, -1475.12481448, -1261.04272428,
        -189.53611144,    24.54597876,   238.62806895,   452.71015915,
        1524.21677199,  1738.29886219,  1952.38095238],[-1702.02087747,   383.86841795,  -545.84986222, -1475.56814238,
        1540.0394332 ,   610.32115304,  -319.39712713, -1249.1154073 ,
        1766.49216829,   836.77388812,   -92.94439205, -1022.66267221,
       -1952.38095238,  1992.94490337,  1063.22662321]])

    
    ind = 3
    coe, _, _ = coe_objective(turbine_x_array[ind], turbine_y_array[ind], solar[ind], battery[3])
    print(coe)
    # coe_array = np.zeros(len(battery))
    # cost_array = np.zeros(len(battery))
    # energy_array = np.zeros(len(battery))

    # for i in range(len(battery)):
    #     print(i)
    #     coe_array[i], cost_array[i], energy_array[i] = coe_objective(turbine_x_array[6], turbine_y_array[6], solar[6], battery[i])
    # print(coe_array)
    # plt.figure(1)
    # plt.plot(coe_array)
    # plt.title("COE")

    # plt.figure(2)
    # plt.plot(cost_array)
    # plt.title("cost")

    # plt.figure(3)
    # plt.plot(energy_array)
    # plt.title("energy")

    # plt.show()
    # print("COE: ", coe)

    # wind_power = hybrid_plant.wind_power
    # solar_power = hybrid_plant.solar_power_with_losses
    # hybrid_power_no_battery = hybrid_plant.hybrid_power_no_battery
    # hybrid_power_with_battery = hybrid_plant.hybrid_power_with_battery

    # plt.figure(1)
    # plt.plot(wind_power, label="wind")
    # plt.plot(solar_power,label="solar")
    # plt.plot(hybrid_power_no_battery,label="wind+solar")
    # plt.plot(hybrid_power_with_battery,label="wind+solar+battery")
    # plt.legend()


    # plt.figure(2)
    # solar_poly = hybrid_plant.solar.solar_geometry
    # plot_poly(solar_poly,ax=plt.gca())
    # plot_turbines(turbine_x,turbine_y,rotor_diameter/2,ax=plt.gca())

    # plt.axis("equal")
    # plt.show()