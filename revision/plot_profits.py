import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interpolate
from plotting_functions import plot_poly, plot_turbines
from make_solar import PlaceSolar
import os
from dotenv import load_dotenv
from hybrid.sites import SiteInfo
from hybrid.sites import flatirons_site as sample_site
from hybrid.keys import set_developer_nrel_gov_key
import json
from shapely.geometry import Polygon
import scipy.interpolate
from simple_dispatch import SimpleDispatch
from hybrid_plant import HybridPlant

def objectives(turbine_x, turbine_y, target_solar_capacity, battery_capacity):
    global hybrid_plant
    global wind_function
    global solar_function
    global battery_function
    global turbine_rating
    global wind_om
    global solar_om
    global battery_om
    global ppa_array
    global interconnect

    if target_solar_capacity == 0.0 and len(turbine_x) == 0:
        return 1E20

    hybrid_plant.update_turbine_locs(turbine_x, turbine_y)
    hybrid_plant.solar_capacity = target_solar_capacity+1E-6
    hybrid_plant.battery_size = battery_capacity

    hybrid_plant.evaluate_plant()
    hybrid_power = hybrid_plant.hybrid_power_with_battery
    print(np.min(hybrid_power))
    # print("wind_power: ", np.sum(hybrid_plant.wind_power))
    # print("solar_power: ", np.sum(hybrid_plant.solar_power_with_losses))
    # yearly_energy = np.sum(hybrid_power)/1E3

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

    # ppa = 70.0
    money_in = hybrid_power*ppa_array/1000.0
    annual_profit = np.sum(money_in) - yearly_cost

    coe = yearly_cost/sum(hybrid_power) * 1000.0
    print("COE: ", coe)
    print("profit: ", annual_profit)
    capacity_factor = np.sum(hybrid_power)/(interconnect*len(hybrid_power))
    print("capacity_factor =", capacity_factor)



if __name__=="__main__":
  rotor_diameter = 200.0
  eval = True

  if eval:
    global hybrid_plant
    global wind_function
    global solar_function
    global battery_function
    global turbine_rating
    global wind_om
    global solar_om
    global battery_om
    global ppa_array
    global interconnect

    ntime = 8760
    # ppa_array = np.zeros(ntime) + 70.0
    time_array = np.arange(ntime)
    ppa_array = np.zeros(ntime)
    for k in range(ntime):
        hour = time_array[k]%24
        # if hour >= 6 and hour < 18:
        #     ppa_array[k] = 70.0
        # else:
        #     ppa_array[k] = 35.0
        if hour < 6 or hour >= 18:
            ppa_array[k] = 70.0
        else:
            ppa_array[k] = 35.0

    # SETUP PLANT EVALUATIONS
    # Set API key
    load_dotenv()
    NREL_API_KEY = os.getenv("NREL_API_KEY")
    set_developer_nrel_gov_key(NREL_API_KEY)  # Set this key manually here if you are not setting it using the .env
      
    outage_start = 3008
    outage_duration = 12
    interconnect = 300000.0
    powercurve_filename = 'turbine_data/high_7r_200d_135h.txt'
    
    hub_height = 135.0
    turbine_rating = 7000.0
    lat = 40.966
    lon = -84.598
    year = 2012
    sample_site['year'] = year
    sample_site['lat'] = lat
    sample_site['lon'] = lon
    site = SiteInfo(sample_site, hub_height=hub_height)
    powercurve_file = open(powercurve_filename)
    powercurve_data = json.load(powercurve_file)
    powercurve_file.close()
    min_turbine_spacing = 4*rotor_diameter
    # corner = min_turbine_spacing*(N-1)/2
    corner = 2500.0
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
    solar_kw_per_km2 = 1000.0/5.0 * 247.105
    # battery
    # min_power = interconnect
    min_power = 10000.0
    battery_dispatch = SimpleDispatch(ntime, min_power)
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


  corner = 2500.0
  ngrid = 25 # dunno what this should be, affects how coarse the solar grid is
  x = np.linspace(-corner,corner,ngrid)
  y = np.linspace(-corner,corner,ngrid)
  xlocs = [[i for i in x] for j in y]
  ylocs = [[j for i in x] for j in x]
  grid_locs = np.zeros((np.shape(xlocs)[0],np.shape(xlocs)[1],2))
  grid_locs[:,:,0] = xlocs[:]
  grid_locs[:,:,1] = ylocs[:]
  grid_locs = np.ndarray.tolist(grid_locs)

  dx = grid_locs[0][1][0] - grid_locs[0][0][0]
  dy = grid_locs[1][0][1] - grid_locs[0][0][1]
  solar_cell_area = dx*dy
  solar_kw_per_km2 = 1000.0/5.0 * 247.105

  min_solar_spacing = 2*rotor_diameter
  place_solar = PlaceSolar(grid_locs,min_solar_spacing)


  plt.figure(figsize=(5.5,6))

  ax1 = plt.subplot(331)
  ax2 = plt.subplot(332)
  ax3 = plt.subplot(333)
  ax4 = plt.subplot(334)
  ax5 = plt.subplot(335)
  ax6 = plt.subplot(336)
  ax7 = plt.subplot(337)
  ax8 = plt.subplot(338)
  ax9 = plt.subplot(339)

  ntime = 1000
  time_arr = np.linspace(0,24,ntime)

  ppa_arr = np.zeros(ntime) + 70.0
  ax1.plot(time_arr, ppa_arr)
  ax1.set_ylim(30,75)
  ax1.set_xlabel("time of day (hr)", fontsize=8)
  ax1.set_ylabel("PPA ($/MWh)",fontsize=8)
  ax1.set_title("constant PPA", fontsize=8)
  ax1.tick_params(axis='both', which='major', labelsize=8)
  ax1.tick_params(axis='both', which='minor', labelsize=8)


  ppa_arr = np.zeros(ntime)
  for k in range(ntime):
      hour = time_arr[k]%24
      if hour >= 6 and hour < 18:
          ppa_arr[k] = 70.0
      else:
          ppa_arr[k] = 35.0

  ax2.plot(time_arr, ppa_arr)
  ax2.set_ylim(30,75)
  ax2.set_xlabel("time of day (hr)", fontsize=8)
  ax2.set_ylabel("PPA ($/MWh)",fontsize=8)
  ax2.set_title("daytime peak", fontsize=8)
  ax2.tick_params(axis='both', which='major', labelsize=8)
  ax2.tick_params(axis='both', which='minor', labelsize=8)



  ppa_arr = np.zeros(ntime)
  for k in range(ntime):
      hour = time_arr[k]%24
      if hour < 6 or hour >= 18:
          ppa_arr[k] = 70.0
      else:
          ppa_arr[k] = 35.0

  ax3.plot(time_arr, ppa_arr)
  ax3.set_ylim(30,75)
  ax3.set_xlabel("time of day (hr)", fontsize=8)
  ax3.set_ylabel("PPA ($/MWh)",fontsize=8)
  ax3.set_title("nighttime peak", fontsize=8)
  ax3.tick_params(axis='both', which='major', labelsize=8)
  ax3.tick_params(axis='both', which='minor', labelsize=8)


  corner = 2600.0
  boundary_polygon = Polygon(((-corner,-corner),(-corner,corner),(corner,corner),(corner,-corner)))
  plot_poly(boundary_polygon,ax=ax4,color="black",fill=False)
  plot_poly(boundary_polygon,ax=ax5,color="black",fill=False)
  plot_poly(boundary_polygon,ax=ax6,color="black",fill=False)
  plot_poly(boundary_polygon,ax=ax7,color="black",fill=False)
  plot_poly(boundary_polygon,ax=ax8,color="black",fill=False)
  plot_poly(boundary_polygon,ax=ax9,color="black",fill=False)

  # # not_resil_const
  # profit = 37077516.13549378
  # wind = 301000.0
  # solar = 42900.173611111626
  # battery = 0.0

  # coe = 36.92644114022708
  # profit = 37077516.39740612

  # turbine_x = np.array([ 2427.92279693,  2105.47794305,  1599.68750281,  2288.82352941,
  #       1783.03308917,  1277.24264893,   771.45220869,  2472.16911577,
  #       1966.37867553,  1460.58823529,   954.79779506,   449.00735482,
  #        -56.78308542,  2149.72426189,  1643.93382165,  1138.14338142,
  #        632.35294118,   126.56250094,  -379.2279393 ,  -885.01837954,
  #       1321.48896778,   815.69852754,   309.9080873 ,  -195.88235294,
  #       -701.67279318, -1207.46323342, -1713.25367366,   493.25367366,
  #        -12.53676658,  -518.32720682, -1024.11764706, -1529.9080873 ,
  #      -2035.69852754,  -334.98162046,  -840.7720607 , -1346.56250094,
  #      -1852.35294118, -2358.14338142, -1163.21691458, -1669.00735482,
  #      -2174.79779506, -1991.45220869, -2497.24264893])

  # turbine_y = np.array([-2470.        , -1641.76470588, -2470.        ,    14.70588235,
  #       -813.52941176, -1641.76470588, -2470.        ,  1671.17647059,
  #        842.94117647,    14.70588235,  -813.52941176, -1641.76470588,
  #      -2470.        ,  2499.41176471,  1671.17647059,   842.94117647,
  #         14.70588235,  -813.52941176, -1641.76470588, -2470.        ,
  #       2499.41176471,  1671.17647059,   842.94117647,    14.70588235,
  #       -813.52941176, -1641.76470588, -2470.        ,  2499.41176471,
  #       1671.17647059,   842.94117647,    14.70588235,  -813.52941176,
  #      -1641.76470588,  2499.41176471,  1671.17647059,   842.94117647,
  #         14.70588235,  -813.52941176,  2499.41176471,  1671.17647059,
  #        842.94117647,  2499.41176471,  1671.17647059])


  profit = 5645419.810419798
  wind = 231000.0
  solar = 296011.1979166696
  battery = 222873.90029325514

  coe = 63.86832137209415
  profit = 5645420.447551072

  turbine_x = np.array([-2326.83074483, -2326.83074483, -2326.83074483, -2326.83074483,
      -2326.83074483, -2326.83074483, -1132.35294118, -1132.35294118,
      -1132.35294118, -1132.35294118, -1132.35294118, -1132.35294118,
          62.12486248,    62.12486248,    62.12486248,    62.12486248,
          62.12486248,    62.12486248,    62.12486248,  1256.60266613,
          1256.60266613,  1256.60266613,  1256.60266613,  1256.60266613,
          1256.60266613,  1256.60266613,  2451.08046978,  2451.08046978,
          2451.08046978,  2451.08046978,  2451.08046978,  2451.08046978,
          2451.08046978])

  turbine_y = np.array([-1761.06069328,  -960.51396801,  -159.96724274,   640.57948252,
          1441.12620779,  2241.67293306, -1716.89278342,  -916.34605815,
          -115.79933289,   684.74739238,  1485.29411765,  2285.84084291,
      -2473.27159883, -1672.72487356,  -872.17814829,   -71.63142303,
          728.91530224,  1529.46202751,  2330.00875277, -2429.10368897,
      -1628.5569637 ,  -828.01023843,   -27.46351317,   773.0832121 ,
          1573.62993736,  2374.17666263, -2384.93577911, -1584.38905384,
          -783.84232858,    16.70439669,   817.25112196,  1617.79784722,
          2418.34457249])
  objectives(turbine_x, turbine_y, solar, battery)
          
  # nturbs = len(turbine_x)
  # turbine_locs = np.zeros((nturbs,2))
  # turbine_locs[:,0] = turbine_x
  # turbine_locs[:,1] = turbine_y
  # place_solar.set_turbine_locs(turbine_locs)
  # place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
  # place_solar.place_solar()
  # solar_poly = place_solar.solar_geometry

  # plot_poly(solar_poly,ax=ax4,color="C1")
  # plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax4)

  # ax4.axis("off")
  # ax4.axis("square")

  # DY = 3200
  # ax4.text(0,0-DY,"COE: %s"%round(coe,2),fontsize=8,horizontalalignment="center")
  # ax4.text(0,-500-DY,"profit: %s"%round(profit/1E6,2),fontsize=8,horizontalalignment="center")
  # ax4.text(0,-1000-DY,"solar: %s"%round(solar/1E3,1),fontsize=8,horizontalalignment="center")
  # ax4.text(0,-1500-DY,"wind: %s"%round(wind/1E3,1),fontsize=8,horizontalalignment="center")
  # ax4.text(0,-2000-DY,"battery: %s"%round(battery/1E3,1),fontsize=8,horizontalalignment="center")

  # # # NOT RESIL DAY
  # profit = 13776322.4929391
  # wind = 301000.0
  # solar = 51480.208333333896
  # battery = 0.0

  # coe = 37.00385988634111
  # profit = 13776323.06809669

  # turbine_x = np.array([-2.49042227e+03, -1.66048144e+03, -2.48532039e+03, -5.70165851e+00,
  #      -8.30540607e+02, -1.65537956e+03, -2.48021850e+03,  8.24239175e+02,
  #      -5.99773294e-01, -8.25438722e+02, -1.65027767e+03, -2.47511662e+03,
  #       2.47901896e+03,  1.65418001e+03,  8.29341060e+02,  4.50211192e+00,
  #      -8.20336837e+02, -1.64517579e+03, -2.47001473e+03,  2.48412084e+03,
  #       1.65928189e+03,  8.34442946e+02,  9.60399714e+00, -8.15234951e+02,
  #      -1.64007390e+03, -2.46491285e+03,  2.48922273e+03,  1.66438378e+03,
  #       8.39544831e+02,  1.47058824e+01, -8.10133066e+02, -1.63497201e+03,
  #      -2.45981096e+03,  2.49432461e+03,  1.66948566e+03,  8.44646716e+02,
  #       1.98077676e+01, -8.05031181e+02,  2.49942650e+03,  1.67458755e+03,
  #       8.49748601e+02,  2.49096528e+01,  1.67968944e+03])
  
  # turbine_y = np.array([-2469.38194084, -2195.06244282, -1641.16236054, -2474.64302708,
  #      -1920.7429448 , -1366.84286252,  -812.94278024, -2200.32352906,
  #      -1646.42344678, -1092.5233645 ,  -538.62328222,    15.27680006,
  #      -2479.90411331, -1926.00403104, -1372.10394876,  -818.20386648,
  #       -264.3037842 ,   289.59629808,   843.49638036, -1651.68453302,
  #      -1097.78445074,  -543.88436846,    10.01571382,   563.9157961 ,
  #       1117.81587838,  1671.71596066,  -823.46495272,  -269.56487044,
  #        284.33521184,   838.23529412,  1392.1353764 ,  1946.03545867,
  #       2499.93554095,     4.75462758,   558.65470986,  1112.55479214,
  #       1666.45487442,  2220.35495669,   832.97420788,  1386.87429016,
  #       1940.77437244,  2494.67445471,  2215.09387046])

  # # objectives(turbine_x, turbine_y, solar, battery)


  # nturbs = len(turbine_x)
  # turbine_locs = np.zeros((nturbs,2))
  # turbine_locs[:,0] = turbine_x
  # turbine_locs[:,1] = turbine_y
  # place_solar.set_turbine_locs(turbine_locs)
  # place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
  # place_solar.place_solar()
  # solar_poly = place_solar.solar_geometry

  # plot_poly(solar_poly,ax=ax5,color="C1")
  # plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax5)

  # ax5.axis("off")
  # ax5.axis("square")

  # DY = 3200
  # ax5.text(0,0-DY,"COE: %s"%round(coe,2),fontsize=8,horizontalalignment="center")
  # ax5.text(0,-500-DY,"profit: %s"%round(profit/1E6,2),fontsize=8,horizontalalignment="center")
  # ax5.text(0,-1000-DY,"solar: %s"%round(solar/1E3,1),fontsize=8,horizontalalignment="center")
  # ax5.text(0,-1500-DY,"wind: %s"%round(wind/1E3,1),fontsize=8,horizontalalignment="center")
  # ax5.text(0,-2000-DY,"battery: %s"%round(battery/1E3,1),fontsize=8,horizontalalignment="center")


  # # NOT RESIL NIGHT
  # profit = 21421051.799398437
  # wind = 301000.0
  # solar = 0.0
  # battery = 0.0

  # coe = 36.79102324715748
  # profit = 21421051.79939837

  # turbine_x = np.array([-1882.51566565,  -906.67110902,    69.1734476 ,  1045.01800423,
  #       2020.86256086, -2376.51607736, -1400.67152073,  -424.82696411,
  #        551.01759252,  1526.86214915, -1894.67193245,  -918.82737582,
  #         57.01718081,  1032.86173744,  2008.70629406, -2388.67234416,
  #      -1412.82778753,  -436.9832309 ,   538.86132573,  1514.70588235,
  #       2490.55043898, -1906.82819924,  -930.98364261,    44.86091402,
  #       1020.70547064,  1996.55002727, -2400.82861095, -1424.98405432,
  #       -449.13949769,   526.70505893,  1502.54961556,  2478.39417219,
  #      -1918.98446603,  -943.13990941,    32.70464722,  1008.54920385,
  #       1984.39376048, -2412.98487774, -1437.14032112,  -461.29576449,
  #        514.54879214,  1490.39334877,  2466.23790539])
  # turbine_y = np.array([ 2486.50177944,  2474.47879432,  2462.4558092 ,  2450.43282408,
  #       2438.40983896,  1790.45125953,  1778.42827441,  1766.40528929,
  #       1754.38230417,  1742.35931905,  1082.3777545 ,  1070.35476938,
  #       1058.33178426,  1046.30879914,  1034.28581402,   386.3272346 ,
  #        374.30424948,   362.28126436,   350.25827924,   338.23529412,
  #        326.212309  ,  -321.74627043,  -333.76925555,  -345.79224067,
  #       -357.81522579,  -369.83821091, -1017.79679034, -1029.81977546,
  #      -1041.84276058, -1053.8657457 , -1065.88873082, -1077.91171594,
  #      -1725.87029536, -1737.89328048, -1749.9162656 , -1761.93925072,
  #      -1773.96223584, -2421.92081527, -2433.94380039, -2445.96678551,
  #      -2457.98977063, -2470.01275575, -2482.03574087])

  # # objectives(turbine_x, turbine_y, solar, battery)

  # nturbs = len(turbine_x)
  # turbine_locs = np.zeros((nturbs,2))
  # turbine_locs[:,0] = turbine_x
  # turbine_locs[:,1] = turbine_y
  # place_solar.set_turbine_locs(turbine_locs)
  # place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
  # place_solar.place_solar()
  # solar_poly = place_solar.solar_geometry

  # plot_poly(solar_poly,ax=ax6,color="C1")
  # plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax6)

  # ax6.axis("off")
  # ax6.axis("square")

  # DY = 3200
  # ax6.text(0,0-DY,"COE: %s"%round(coe,2),fontsize=8,horizontalalignment="center")
  # ax6.text(0,-500-DY,"profit: %s"%round(profit/1E6,2),fontsize=8,horizontalalignment="center")
  # ax6.text(0,-1000-DY,"solar: %s"%round(solar/1E3,1),fontsize=8,horizontalalignment="center")
  # ax6.text(0,-1500-DY,"wind: %s"%round(wind/1E3,1),fontsize=8,horizontalalignment="center")
  # ax6.text(0,-2000-DY,"battery: %s"%round(battery/1E3,1),fontsize=8,horizontalalignment="center")

  # # # RESIL CONST
  # profit = 5645419.810419798
  # wind = 231000.0
  # solar = 296011.1979166696
  # battery = 222873.90029325514

  # coe = 63.86832137209415
  # profit = 5645420.447551072

  # turbine_x = np.array([-2326.83074483, -2326.83074483, -2326.83074483, -2326.83074483,
  #      -2326.83074483, -2326.83074483, -1132.35294118, -1132.35294118,
  #      -1132.35294118, -1132.35294118, -1132.35294118, -1132.35294118,
  #         62.12486248,    62.12486248,    62.12486248,    62.12486248,
  #         62.12486248,    62.12486248,    62.12486248,  1256.60266613,
  #       1256.60266613,  1256.60266613,  1256.60266613,  1256.60266613,
  #       1256.60266613,  1256.60266613,  2451.08046978,  2451.08046978,
  #       2451.08046978,  2451.08046978,  2451.08046978,  2451.08046978,
  #       2451.08046978])

  # turbine_y = np.array([-1761.06069328,  -960.51396801,  -159.96724274,   640.57948252,
  #       1441.12620779,  2241.67293306, -1716.89278342,  -916.34605815,
  #       -115.79933289,   684.74739238,  1485.29411765,  2285.84084291,
  #      -2473.27159883, -1672.72487356,  -872.17814829,   -71.63142303,
  #        728.91530224,  1529.46202751,  2330.00875277, -2429.10368897,
  #      -1628.5569637 ,  -828.01023843,   -27.46351317,   773.0832121 ,
  #       1573.62993736,  2374.17666263, -2384.93577911, -1584.38905384,
  #       -783.84232858,    16.70439669,   817.25112196,  1617.79784722,
  #       2418.34457249])

  # # objectives(turbine_x, turbine_y, solar, battery)
  
  # nturbs = len(turbine_x)
  # turbine_locs = np.zeros((nturbs,2))
  # turbine_locs[:,0] = turbine_x
  # turbine_locs[:,1] = turbine_y
  # place_solar.set_turbine_locs(turbine_locs)
  # place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
  # place_solar.place_solar()
  # solar_poly = place_solar.solar_geometry

  # plot_poly(solar_poly,ax=ax7,color="C1")
  # plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax7)

  # ax7.axis("off")
  # ax7.axis("square")

  # DY = 3200
  # ax7.text(0,0-DY,"COE: %s"%round(coe,2),fontsize=8,horizontalalignment="center")
  # ax7.text(0,-500-DY,"profit: %s"%round(profit/1E6,2),fontsize=8,horizontalalignment="center")
  # ax7.text(0,-1000-DY,"solar: %s"%round(solar/1E3,1),fontsize=8,horizontalalignment="center")
  # ax7.text(0,-1500-DY,"wind: %s"%round(wind/1E3,1),fontsize=8,horizontalalignment="center")
  # ax7.text(0,-2000-DY,"battery: %s"%round(battery/1E3,1),fontsize=8,horizontalalignment="center")
  
  # # # RESIL DAY
  # profit = -11767038.269357637
  # wind = 182000.0
  # solar = 555557.2482638895
  # battery = 197458.45552297164

  # coe = 65.88958604373684
  # profit = -11767029.861832589

  # turbine_x = np.array([-1980.50573774,  -531.12147222,   918.2627933 ,  2367.64705882,
  #      -1975.40385253,  -526.019587  ,   923.36467852,  2372.74894404,
  #      -1970.30196731,  -520.91770179,   928.46656373,  2377.85082925,
  #      -1965.2000821 ,  -515.81581657,   933.56844895,  2382.95271447,
  #      -1960.09819688,  -510.71393136,   938.67033416,  2388.05459969,
  #      -1954.99631166,  -505.61204614,   943.77221938,  2393.1564849 ,
  #        948.8741046 ,  2398.25837012])
  
  # turbine_y = np.array([ 2.45850922e+03,  2.46743752e+03,  2.47636582e+03,  2.48529412e+03,
  #       1.63028964e+03,  1.63921794e+03,  1.64814624e+03,  1.65707454e+03,
  #       8.02070060e+02,  8.10998359e+02,  8.19926658e+02,  8.28854957e+02,
  #      -2.61495206e+01, -1.72212215e+01, -8.29292237e+00,  6.35376752e-01,
  #      -8.54369101e+02, -8.45440802e+02, -8.36512503e+02, -8.27584204e+02,
  #      -1.68258868e+03, -1.67366038e+03, -1.66473208e+03, -1.65580378e+03,
  #      -2.49295166e+03, -2.48402336e+03])
  # # objectives(turbine_x, turbine_y, solar, battery)
  
  # nturbs = len(turbine_x)
  # turbine_locs = np.zeros((nturbs,2))
  # turbine_locs[:,0] = turbine_x
  # turbine_locs[:,1] = turbine_y
  # place_solar.set_turbine_locs(turbine_locs)
  # place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
  # place_solar.place_solar()
  # solar_poly = place_solar.solar_geometry

  # plot_poly(solar_poly,ax=ax8,color="C1")
  # plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax8)

  # ax8.axis("off")
  # ax8.axis("square")

  # DY = 3200
  # ax8.text(0,0-DY,"COE: %s"%round(coe,2),fontsize=8,horizontalalignment="center")
  # ax8.text(0,-500-DY,"profit: %s"%round(profit/1E6,2),fontsize=8,horizontalalignment="center")
  # ax8.text(0,-1000-DY,"solar: %s"%round(solar/1E3,1),fontsize=8,horizontalalignment="center")
  # ax8.text(0,-1500-DY,"wind: %s"%round(wind/1E3,1),fontsize=8,horizontalalignment="center")
  # ax8.text(0,-2000-DY,"battery: %s"%round(battery/1E3,1),fontsize=8,horizontalalignment="center")

  # # # # RESIL NIGHT


  # """from the optimizer for night"""
  # # profit = -8940632.54685422
  # # wind = 203000.0
  # # solar = 360361.4583333376
  # # battery = 216031.28054740958

  # # coe = 65.65618195635238
  # # profit = -8940637.954323769

  # # turbine_x = np.array([ 1921.79854451,  2318.60683387,   485.79518599,   882.60347535,
  # #       1279.41176471,  1676.22005406,  2073.02834342,  2469.83663278,
  # #       -950.20817253,  -553.39988317,  -156.59159381,   240.21669554,
  # #        637.0249849 ,  1033.83327426,  1430.64156362,  1827.44985297,
  # #      -2386.21153105, -1989.40324169, -1592.59495233, -1195.78666298,
  # #       -798.97837362,  -402.17008426,    -5.3617949 ,   391.44649445,
  # #      -2234.98173214, -1838.17344278, -1441.36515342, -1044.55686407,
  # #      -2480.56022259])

  # # turbine_y = np.array([ 2497.96180169,  1800.78678054,  2497.29121876,  1800.11619762,
  # #       1102.94117647,   405.76615532,  -291.40886582,  -988.58388697,
  # #       2496.62063584,  1799.44561469,  1102.27059355,   405.0955724 ,
  # #       -292.07944875,  -989.25446989, -1686.42949104, -2383.60451219,
  # #       2495.95005292,  1798.77503177,  1101.60001062,   404.42498947,
  # #       -292.75003167,  -989.92505282, -1687.10007397, -2384.27509511,
  # #       -293.4206146 ,  -990.59563574, -1687.77065689, -2384.94567804,
  # #      -2385.61626096])
  # # objectives(turbine_x, turbine_y, solar, battery)

  # """from the optimizer for const (which is better for this too)"""
  # wind = 231000.0
  # solar = 296011.1979166696
  # battery = 222873.90029325514

  # coe = 63.86832137209415
  # profit = -8054800.

  # turbine_x = np.array([-2326.83074483, -2326.83074483, -2326.83074483, -2326.83074483,
  #      -2326.83074483, -2326.83074483, -1132.35294118, -1132.35294118,
  #      -1132.35294118, -1132.35294118, -1132.35294118, -1132.35294118,
  #         62.12486248,    62.12486248,    62.12486248,    62.12486248,
  #         62.12486248,    62.12486248,    62.12486248,  1256.60266613,
  #       1256.60266613,  1256.60266613,  1256.60266613,  1256.60266613,
  #       1256.60266613,  1256.60266613,  2451.08046978,  2451.08046978,
  #       2451.08046978,  2451.08046978,  2451.08046978,  2451.08046978,
  #       2451.08046978])

  # turbine_y = np.array([-1761.06069328,  -960.51396801,  -159.96724274,   640.57948252,
  #       1441.12620779,  2241.67293306, -1716.89278342,  -916.34605815,
  #       -115.79933289,   684.74739238,  1485.29411765,  2285.84084291,
  #      -2473.27159883, -1672.72487356,  -872.17814829,   -71.63142303,
  #        728.91530224,  1529.46202751,  2330.00875277, -2429.10368897,
  #      -1628.5569637 ,  -828.01023843,   -27.46351317,   773.0832121 ,
  #       1573.62993736,  2374.17666263, -2384.93577911, -1584.38905384,
  #       -783.84232858,    16.70439669,   817.25112196,  1617.79784722,
  #       2418.34457249])

  # nturbs = len(turbine_x)
  # turbine_locs = np.zeros((nturbs,2))
  # turbine_locs[:,0] = turbine_x
  # turbine_locs[:,1] = turbine_y
  # place_solar.set_turbine_locs(turbine_locs)
  # place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
  # place_solar.place_solar()
  # solar_poly = place_solar.solar_geometry

  # plot_poly(solar_poly,ax=ax9,color="C1")
  # plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax9)

  # ax9.axis("off")
  # ax9.axis("square")

  # DY = 3200
  # ax9.text(0,0-DY,"COE: %s"%round(coe,2),fontsize=8,horizontalalignment="center")
  # ax9.text(0,-500-DY,"profit: %s"%round(profit/1E6,2),fontsize=8,horizontalalignment="center")
  # ax9.text(0,-1000-DY,"solar: %s"%round(solar/1E3,1),fontsize=8,horizontalalignment="center")
  # ax9.text(0,-1500-DY,"wind: %s"%round(wind/1E3,1),fontsize=8,horizontalalignment="center")
  # ax9.text(0,-2000-DY,"battery: %s"%round(battery/1E3,1),fontsize=8,horizontalalignment="center")


  # ax4.text(-3500,0,"no minimum\npower constraint",fontsize=8,rotation=90,
  #     horizontalalignment="center",verticalalignment="center")
  # ax7.text(-3500,0,"10-MW minimum\npower constraint",fontsize=8,rotation=90,
  #     horizontalalignment="center",verticalalignment="center")

  # plt.subplots_adjust(left=0.08,right=0.99,top=0.95,bottom=0.1,wspace=0.4,hspace=0.45)
  
  # plt.savefig("make_figures/figures/ppa_layouts_revision.pdf",transparent=True)
  # plt.show()