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
    outage_duration = 30
    interconnect = 100000.0
    min_power = 10000.0
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

    
    solar = np.array([223421.25207756367,201517.2077562349,365797.5401662054,357035.9224376759,
                  118281.83933518121,374559.15789473767,359226.3268698072,
                  127043.45706371329,129233.86149584601,181803.56786703764,
                  359226.32686980703,249706.10526316008])
    battery = np.array([250244.37927663734,250244.37927663734,250244.37927663734,250244.37927663734,
                    250244.37927663734,250244.37927663734,250244.37927663734,
                    250244.37927663734,250244.37927663734,265884.65298142715,
                    281524.926686217,312805.4740957967])

    turbine_x_array = np.array([[ 1905.87477526,  1809.52380952,  1098.57670674,  1002.225741  ,
         905.87477526,   809.52380952,    98.57670674,     2.225741  ,
         -94.12522474,  -190.47619048,  -901.42329326,  -997.774259  ,
       -1094.12522474, -1190.47619048, -1901.42329326, -1997.774259  ],[ 1913.98542614,  1845.98766558,  1777.98990502,   997.0689027 ,
         929.07114213,   861.07338157,    12.15461869,   -55.84314188,
        -123.84090244,  -904.76190476,  -972.75966532, -1040.75742589,
       -1821.67842821, -1889.67618877, -1957.67394933],[-1453.44340311,  -411.60773749,   630.22792812,  1672.06359374,
       -1043.96056707,    -2.12490145,  1039.71076416, -1676.31339665,
        -634.47773103,   407.35793458,  1449.1936002 , -1266.83056061,
        -224.99489499,   816.84077062,  1858.67643624],[-1890.2437802 , -1029.07161889, -1913.64237371, -1052.4702124 ,
        -191.29805109,   669.87411022, -1937.04096721, -1075.8688059 ,
        -214.69664459,   646.47551672,  1507.64767803,  -238.0952381 ,
         623.07692322,  1484.24908453,  1460.85049102],[ 1952.38095238,  1923.8798746 ,  1895.37879683,  1866.87771905,
         666.66666667,   638.16558889,   609.66451111,   581.16343334,
        -619.04761905,  -647.54869682,  -676.0497746 ,  -704.55085238,
       -1904.76190476, -1933.26298254, -1961.76406032, -1990.26513809],[-1456.5340372 ,    67.72015634,  1591.97434987,  -605.79510314,
         918.45909039, -1279.31036262,   244.94383091,  1769.19802445,
       -1952.82562211,  -428.57142857,  1095.68276496, -1102.08668805,
         422.16750548,  1946.42169902],[-1953.98043925,  -667.4664101 ,   619.04761905,  1905.56164819,
       -1925.48822012,  -638.97419097,   647.53983818,  1934.05386732,
       -1896.99600099,  -610.48197184,   676.03205731,  1962.54608645,
       -1868.50378186,  -581.98975271,   704.52427644,  1991.03830558],[-1872.17753173,  -603.34163096, -1038.25799997, -1473.17436897,
       -1908.09073798,   665.49426981,   230.5779008 ,  -204.3384682 ,
        -639.25483721,  1934.33017057,  1499.41380157,  1064.49743256,
         629.58106356,  1898.41696433],[ 1729.44325771,  1801.89995477,  1874.35665183,  1946.81334889,
         549.10215154,   621.5588486 ,   694.01554566,  -703.69565169,
        -631.23895463,  -558.78225756, -1956.49345492, -1884.03675785,
       -1811.58006079, -1739.12336373],[-1958.63614107,  -797.2276643 ,   364.18081247, -1488.59373749,
        -327.18526072,   834.22321605,  1995.63169282, -1018.55133391,
         142.85714286,  1304.26561963, -1709.9174071 ,  -548.50893033,
         612.89954643,  1774.3080232 ],[ 1952.38095238,  1952.38095238,  1952.38095238,  1952.38095238,
         666.66666667,   666.66666667,   666.66666667,   666.66666667,
        -619.04761905,  -619.04761905,  -619.04761905,  -619.04761905,
       -1904.76190476, -1904.76190476, -1904.76190476, -1904.76190476],[ 1667.11706571,   723.67916702,  1735.07560084,   791.63770215,
        -151.80019654, -1095.23809524,  1803.03413597,   859.59623728,
         -83.84166141, -1027.27956011, -1970.7174588 ,   927.55477241,
         -15.88312629,  -959.32102498, -1902.75892367,  -891.36248985,
       -1834.80038854]])

    turbine_y_array = np.array([[ -666.66666667, -1952.38095238,  1904.76190476,   619.04761905,
        -666.66666667, -1952.38095238,  1904.76190476,   619.04761905,
        -666.66666667, -1952.38095238,  1904.76190476,   619.04761905,
        -666.66666667, -1952.38095238,  1904.76190476,   619.04761905],[ 1135.34764803,  -227.11653711, -1589.58072225,  1671.37062072,
         308.90643558, -1053.55774956,   844.92940827,  -517.53477688,
       -1879.99896202,  1380.95238095,    18.48819581, -1343.97598933,
        1916.97535364,   554.5111685 ,  -807.95301664],[ 1820.21900943,  1859.20174965,  1898.18448988,  1937.1672301 ,
         601.45527865,   640.43801887,   679.42075909,  -656.29119236,
        -617.30845214,  -578.32571192,  -539.3429717 , -1875.05492315,
       -1836.07218292, -1797.0894427 , -1758.10670248],[ 1139.62880443,  1899.61774602,  -397.7004089 ,   362.28853269,
        1122.27747429,  1882.26641588, -1935.02962224, -1175.04068064,
        -415.05173905,   344.93720255,  1104.92614414, -1952.38095238,
       -1192.39201079,  -432.40306919, -1969.73228252],[ 1761.9047619 ,   619.04761905,  -523.80952381, -1666.66666667,
        1761.9047619 ,   619.04761905,  -523.80952381, -1666.66666667,
        1761.9047619 ,   619.04761905,  -523.80952381, -1666.66666667,
        1761.9047619 ,   619.04761905,  -523.80952381, -1666.66666667],[ 1987.25130669,  1987.25130669,  1987.25130669,  1023.24690287,
        1023.24690287,    59.24249905,    59.24249905,    59.24249905,
        -904.76190476,  -904.76190476,  -904.76190476, -1868.76630858,
       -1868.76630858, -1868.76630858],[ 1508.26914618,  1492.22981119,  1476.19047619,  1460.15114119,
         365.7672239 ,   349.7278889 ,   333.6885539 ,   317.64921891,
        -776.73469839,  -792.77403339,  -808.81336838,  -824.85270338,
       -1919.23662068, -1935.27595567, -1951.31529067, -1967.35462567],[-1994.7766471 , -1943.014046  ,  -834.86521037,   273.28362527,
        1381.43246091, -1891.2514449 ,  -783.10260927,   325.04622637,
        1433.19506201, -1839.4888438 ,  -731.34000817,   376.80882747,
        1484.95766311,  1536.72026421],[ 1926.51250519,   765.60875797,  -395.29498925, -1556.19873647,
        1121.7695651 ,   -39.13418212, -1200.03792934,  1477.93037223,
         317.02662501,  -843.87712221,  1834.09117936,   673.18743214,
        -487.71631508, -1648.6200623 ],[ 1450.08286949,  1654.87051957,  1859.65816965,   408.36189542,
         613.1495455 ,   817.93719558,  1022.72484565,  -633.35907865,
        -428.57142857,  -223.78377849, -1879.8677028 , -1675.08005272,
       -1470.29240264, -1265.50475256],[ 1476.19047619,   333.33333333,  -809.52380952, -1952.38095238,
        1476.19047619,   333.33333333,  -809.52380952, -1952.38095238,
        1476.19047619,   333.33333333,  -809.52380952, -1952.38095238,
        1476.19047619,   333.33333333,  -809.52380952, -1952.38095238],[-1.46898138e+03, -1.98275028e+03, -2.20598072e+02, -7.34366969e+02,
       -1.24813587e+03, -1.76190476e+03,  1.02778524e+03,  5.14016339e+02,
        2.47442165e-01, -5.13521454e+02, -1.02729035e+03,  1.76239965e+03,
        1.24863075e+03,  7.34861853e+02,  2.21092957e+02,  1.98324516e+03,
        1.46947626e+03]])

    
    coe, _, _ = coe_objective(turbine_x_array[6], turbine_y_array[6], solar[6], battery[-1])
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