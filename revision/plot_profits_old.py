import numpy as np
import matplotlib.pyplot as plt
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

    ntime = 8760
    ppa_array = np.zeros(ntime) + 70.0
    time_array = np.arange(ntime)
    # ppa_array = np.zeros(ntime)
    # for k in range(ntime):
    #     hour = time_array[k]%24
    #     # if hour >= 6 and hour < 18:
    #     #     ppa_array[k] = 70.0
    #     # else:
    #     #     ppa_array[k] = 35.0
    #     if hour < 6 or hour >= 18:
    #         ppa_array[k] = 70.0
    #     else:
            # ppa_array[k] = 35.0

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
    min_power = interconnect
    # min_power = 10000.0
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


  plt.figure(figsize=(6,6))

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





  # not_resil_const
  profit = 86838641397.65564
  wind = 343000.0
  solar = 75075.30381944594
  battery = 1000000.0
  coe = 125.13824332813539
  profit = 86838641288.36745

  turbine_x = np.array([ 2499.41176471,  2499.41176471,  1671.17647059,  2499.41176471,
          1671.17647059,   842.94117647,  2499.41176471,  1671.17647059,
          842.94117647,    14.70588235,  2499.41176471,  1671.17647059,
          842.94117647,    14.70588235,  -813.52941176,  2499.41176471,
          1671.17647059,   842.94117647,    14.70588235,  -813.52941176,
        -1641.76470588,  2499.41176471,  1671.17647059,   842.94117647,
            14.70588235,  -813.52941176, -1641.76470588, -2470.        ,
          1671.17647059,   842.94117647,    14.70588235,  -813.52941176,
        -1641.76470588, -2470.        ,   842.94117647,    14.70588235,
          -813.52941176, -1641.76470588, -2470.        ,    14.70588235,
          -813.52941176, -1641.76470588, -2470.        ,  -813.52941176,
        -1641.76470588, -2470.        , -1641.76470588, -2470.        ,
        -2470.        ])

  turbine_y = np.array([-2499.41176471, -1671.17647059, -2499.41176471,  -842.94117647,
        -1671.17647059, -2499.41176471,   -14.70588235,  -842.94117647,
        -1671.17647059, -2499.41176471,   813.52941176,   -14.70588235,
          -842.94117647, -1671.17647059, -2499.41176471,  1641.76470588,
          813.52941176,   -14.70588235,  -842.94117647, -1671.17647059,
        -2499.41176471,  2470.        ,  1641.76470588,   813.52941176,
          -14.70588235,  -842.94117647, -1671.17647059, -2499.41176471,
          2470.        ,  1641.76470588,   813.52941176,   -14.70588235,
          -842.94117647, -1671.17647059,  2470.        ,  1641.76470588,
          813.52941176,   -14.70588235,  -842.94117647,  2470.        ,
          1641.76470588,   813.52941176,   -14.70588235,  2470.        ,
          1641.76470588,   813.52941176,  2470.        ,  1641.76470588,
          2470.        ])

  objectives(turbine_x, turbine_y, solar, battery)
          
  nturbs = len(turbine_x)
  turbine_locs = np.zeros((nturbs,2))
  turbine_locs[:,0] = turbine_x
  turbine_locs[:,1] = turbine_y
  place_solar.set_turbine_locs(turbine_locs)
  place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
  place_solar.place_solar()
  solar_poly = place_solar.solar_geometry

  plot_poly(solar_poly,ax=ax4,color="C1")
  plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax4)

  ax4.axis("off")
  ax4.axis("square")



  # NOT RESIL DAY
  profit = 61588078299.62448
  wind = 343000.0
  solar = 77220.31250000153
  battery = 1000000.0

  coe = 125.12885967972228
  profit = 61588079233.07932

  turbine_x = np.array([ 2490.        ,  2490.        ,  2490.        ,  2490.        ,
        2490.        ,  2490.        ,  2490.        ,  1661.76470588,
        1661.76470588,  1661.76470588,  1661.76470588,  1661.76470588,
        1661.76470588,  1661.76470588,   833.52941176,   833.52941176,
         833.52941176,   833.52941176,   833.52941176,   833.52941176,
         833.52941176,     5.29411765,     5.29411765,     5.29411765,
           5.29411765,     5.29411765,     5.29411765,     5.29411765,
        -822.94117647,  -822.94117647,  -822.94117647,  -822.94117647,
        -822.94117647,  -822.94117647,  -822.94117647, -1651.17647059,
       -1651.17647059, -1651.17647059, -1651.17647059, -1651.17647059,
       -1651.17647059, -1651.17647059, -2479.41176471, -2479.41176471,
       -2479.41176471, -2479.41176471, -2479.41176471, -2479.41176471,
       -2479.41176471])
  turbine_y = np.array([ 2494.70588235,  1666.47058824,   838.23529412,    10.        ,
        -818.23529412, -1646.47058824, -2474.70588235,  2494.70588235,
        1666.47058824,   838.23529412,    10.        ,  -818.23529412,
       -1646.47058824, -2474.70588235,  2494.70588235,  1666.47058824,
         838.23529412,    10.        ,  -818.23529412, -1646.47058824,
       -2474.70588235,  2494.70588235,  1666.47058824,   838.23529412,
          10.        ,  -818.23529412, -1646.47058824, -2474.70588235,
        2494.70588235,  1666.47058824,   838.23529412,    10.        ,
        -818.23529412, -1646.47058824, -2474.70588235,  2494.70588235,
        1666.47058824,   838.23529412,    10.        ,  -818.23529412,
       -1646.47058824, -2474.70588235,  2494.70588235,  1666.47058824,
         838.23529412,    10.        ,  -818.23529412, -1646.47058824,
       -2474.70588235])

  # objectives(turbine_x, turbine_y, solar, battery)

  nturbs = len(turbine_x)
  turbine_locs = np.zeros((nturbs,2))
  turbine_locs[:,0] = turbine_x
  turbine_locs[:,1] = turbine_y
  place_solar.set_turbine_locs(turbine_locs)
  place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
  place_solar.place_solar()
  solar_poly = place_solar.solar_geometry

  plot_poly(solar_poly,ax=ax5,color="C1")
  plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax5)

  ax5.axis("off")
  ax5.axis("square")




  # NOT RESIL NIGHT
  profit = 68219479263.066376
  wind = 343000.0
  solar = 27885.11284722276
  battery = 1000000.0

  coe = 125.73682431931252
  profit = 68219478863.8507

  turbine_x = np.array([-2.47022494e+03, -2.47515299e+03, -2.48008104e+03, -2.48500909e+03,
       -2.48993715e+03, -2.49486520e+03, -2.49979325e+03, -1.64205250e+03,
       -1.64698055e+03, -1.65190860e+03, -1.65683666e+03, -1.66176471e+03,
       -1.66669276e+03, -1.67162081e+03, -8.13880063e+02, -8.18808114e+02,
       -8.23736165e+02, -8.28664216e+02, -8.33592266e+02, -8.38520317e+02,
       -8.43448368e+02,  1.42923762e+01,  9.36432542e+00,  4.43627461e+00,
       -4.91776201e-01, -5.41982701e+00, -1.03478778e+01, -1.52759286e+01,
        8.42464816e+02,  8.37536765e+02,  8.32608714e+02,  8.27680663e+02,
        8.22752612e+02,  8.17824562e+02,  8.12896511e+02,  1.67063726e+03,
        1.66570920e+03,  1.66078115e+03,  1.65585310e+03,  1.65092505e+03,
        1.64599700e+03,  1.64106895e+03,  2.49880969e+03,  2.49388164e+03,
        2.48895359e+03,  2.48402554e+03,  2.47909749e+03,  2.47416944e+03,
        2.46924139e+03])
  turbine_y = np.array([-2489.61534154, -1689.61534154,  -889.61534154,   -89.61534154,
         710.38465846,  1510.38465846,  2310.38465846, -2479.41176471,
       -1679.41176471,  -879.41176471,   -79.41176471,   720.58823529,
        1520.58823529,  2320.58823529, -2469.20818787, -1669.20818787,
        -869.20818787,   -69.20818787,   730.79181213,  1530.79181213,
        2330.79181213, -2459.00461103, -1659.00461103,  -859.00461103,
         -59.00461103,   740.99538897,  1540.99538897,  2340.99538897,
       -2448.80103419, -1648.80103419,  -848.80103419,   -48.80103419,
         751.19896581,  1551.19896581,  2351.19896581, -2438.59745735,
       -1638.59745735,  -838.59745735,   -38.59745735,   761.40254265,
        1561.40254265,  2361.40254265, -2428.39388051, -1628.39388051,
        -828.39388051,   -28.39388051,   771.60611949,  1571.60611949,
        2371.60611949])

  # objectives(turbine_x, turbine_y, solar, battery)

  nturbs = len(turbine_x)
  turbine_locs = np.zeros((nturbs,2))
  turbine_locs[:,0] = turbine_x
  turbine_locs[:,1] = turbine_y
  place_solar.set_turbine_locs(turbine_locs)
  place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
  place_solar.place_solar()
  solar_poly = place_solar.solar_geometry

  plot_poly(solar_poly,ax=ax6,color="C1")
  plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax6)

  ax6.axis("off")
  ax6.axis("square")

  # RESIL CONST
  profit = 82403734079.93352
  wind = 343000.0
  solar = 77220.31250000153
  battery = 344086.0215053763

  coe = 72.03971623911045
  profit = 82403733663.15857

  turbine_x = np.array([ 2494.70588235,  2494.70588235,  2494.70588235,  2494.70588235,
        2494.70588235,  2494.70588235,  2494.70588235,  1666.47058824,
        1666.47058824,  1666.47058824,  1666.47058824,  1666.47058824,
        1666.47058824,  1666.47058824,   838.23529412,   838.23529412,
         838.23529412,   838.23529412,   838.23529412,   838.23529412,
         838.23529412,    10.        ,    10.        ,    10.        ,
          10.        ,    10.        ,    10.        ,    10.        ,
        -818.23529412,  -818.23529412,  -818.23529412,  -818.23529412,
        -818.23529412,  -818.23529412,  -818.23529412, -1646.47058824,
       -1646.47058824, -1646.47058824, -1646.47058824, -1646.47058824,
       -1646.47058824, -1646.47058824, -2474.70588235, -2474.70588235,
       -2474.70588235, -2474.70588235, -2474.70588235, -2474.70588235,
       -2474.70588235])
  turbine_y = np.array([ 2.48529412e+03,  1.65705882e+03,  8.28823529e+02,  5.88235294e-01,
       -8.27647059e+02, -1.65588235e+03, -2.48411765e+03,  2.48529412e+03,
        1.65705882e+03,  8.28823529e+02,  5.88235294e-01, -8.27647059e+02,
       -1.65588235e+03, -2.48411765e+03,  2.48529412e+03,  1.65705882e+03,
        8.28823529e+02,  5.88235294e-01, -8.27647059e+02, -1.65588235e+03,
       -2.48411765e+03,  2.48529412e+03,  1.65705882e+03,  8.28823529e+02,
        5.88235294e-01, -8.27647059e+02, -1.65588235e+03, -2.48411765e+03,
        2.48529412e+03,  1.65705882e+03,  8.28823529e+02,  5.88235294e-01,
       -8.27647059e+02, -1.65588235e+03, -2.48411765e+03,  2.48529412e+03,
        1.65705882e+03,  8.28823529e+02,  5.88235294e-01, -8.27647059e+02,
       -1.65588235e+03, -2.48411765e+03,  2.48529412e+03,  1.65705882e+03,
        8.28823529e+02,  5.88235294e-01, -8.27647059e+02, -1.65588235e+03,
       -2.48411765e+03])

  # objectives(turbine_x, turbine_y, solar, battery)
  
  nturbs = len(turbine_x)
  turbine_locs = np.zeros((nturbs,2))
  turbine_locs[:,0] = turbine_x
  turbine_locs[:,1] = turbine_y
  place_solar.set_turbine_locs(turbine_locs)
  place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
  place_solar.place_solar()
  solar_poly = place_solar.solar_geometry

  plot_poly(solar_poly,ax=ax7,color="C1")
  plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax7)

  ax7.axis("off")
  ax7.axis("square")
  
  # RESIL DAY
  profit = 58077934318.80117
  wind = 343000.0
  solar = 77220.31250000153
  battery = 375366.568914956

  coe = 74.85291258412869
  profit = 58077935354.46238

  turbine_x = np.array([ 2494.70588235,  2494.70588235,  2494.70588235,  2494.70588235,
        2494.70588235,  2494.70588235,  2494.70588235,  1666.47058824,
        1666.47058824,  1666.47058824,  1666.47058824,  1666.47058824,
        1666.47058824,  1666.47058824,   838.23529412,   838.23529412,
         838.23529412,   838.23529412,   838.23529412,   838.23529412,
         838.23529412,    10.        ,    10.        ,    10.        ,
          10.        ,    10.        ,    10.        ,    10.        ,
        -818.23529412,  -818.23529412,  -818.23529412,  -818.23529412,
        -818.23529412,  -818.23529412,  -818.23529412, -1646.47058824,
       -1646.47058824, -1646.47058824, -1646.47058824, -1646.47058824,
       -1646.47058824, -1646.47058824, -2474.70588235, -2474.70588235,
       -2474.70588235, -2474.70588235, -2474.70588235, -2474.70588235,
       -2474.70588235])
  
  turbine_y = np.array([ 2.48529412e+03,  1.65705882e+03,  8.28823529e+02,  5.88235294e-01,
       -8.27647059e+02, -1.65588235e+03, -2.48411765e+03,  2.48529412e+03,
        1.65705882e+03,  8.28823529e+02,  5.88235294e-01, -8.27647059e+02,
       -1.65588235e+03, -2.48411765e+03,  2.48529412e+03,  1.65705882e+03,
        8.28823529e+02,  5.88235294e-01, -8.27647059e+02, -1.65588235e+03,
       -2.48411765e+03,  2.48529412e+03,  1.65705882e+03,  8.28823529e+02,
        5.88235294e-01, -8.27647059e+02, -1.65588235e+03, -2.48411765e+03,
        2.48529412e+03,  1.65705882e+03,  8.28823529e+02,  5.88235294e-01,
       -8.27647059e+02, -1.65588235e+03, -2.48411765e+03,  2.48529412e+03,
        1.65705882e+03,  8.28823529e+02,  5.88235294e-01, -8.27647059e+02,
       -1.65588235e+03, -2.48411765e+03,  2.48529412e+03,  1.65705882e+03,
        8.28823529e+02,  5.88235294e-01, -8.27647059e+02, -1.65588235e+03,
       -2.48411765e+03])
  # objectives(turbine_x, turbine_y, solar, battery)
  
  nturbs = len(turbine_x)
  turbine_locs = np.zeros((nturbs,2))
  turbine_locs[:,0] = turbine_x
  turbine_locs[:,1] = turbine_y
  place_solar.set_turbine_locs(turbine_locs)
  place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
  place_solar.place_solar()
  solar_poly = place_solar.solar_geometry

  plot_poly(solar_poly,ax=ax8,color="C1")
  plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax8)

  ax8.axis("off")
  ax8.axis("square")

  # RESIL NIGHT

  profit = 65128391799.89735
  wind = 343000.0
  solar = 25740.104166667184
  battery = 500488.7585532747

  coe = 86.26592773438696
  profit = 65128391230.68032

  turbine_x = np.array([ 2.45983240e+03,  1.65969577e+03,  8.59559151e+02,  5.94225283e+01,
       -7.40714095e+02, -1.54085072e+03, -2.34098734e+03,  2.44452752e+03,
        1.64439089e+03,  8.44254270e+02,  4.41176471e+01, -7.56018976e+02,
       -1.55615560e+03, -2.35629222e+03,  2.42922263e+03,  1.62908601e+03,
        8.28949389e+02,  2.88127658e+01, -7.71323857e+02, -1.57146048e+03,
       -2.37159710e+03,  2.41391775e+03,  1.61378113e+03,  8.13644507e+02,
        1.35078845e+01, -7.86628738e+02, -1.58676536e+03, -2.38690198e+03,
        2.39861287e+03,  1.59847625e+03,  7.98339626e+02, -1.79699679e+00,
       -8.01933620e+02, -1.60207024e+03, -2.40220687e+03,  2.38330799e+03,
        1.58317137e+03,  7.83034745e+02, -1.71018781e+01, -8.17238501e+02,
       -1.61737512e+03, -2.41751175e+03,  2.36800311e+03,  1.56786649e+03,
        7.67729864e+02, -3.24067594e+01, -8.32543382e+02, -1.63268001e+03,
       -2.43281663e+03])

  turbine_y = np.array([-2489.8585792 , -2489.8585792 , -2489.8585792 , -2489.8585792 ,
       -2489.8585792 , -2489.8585792 , -2489.8585792 , -1661.76470588,
       -1661.76470588, -1661.76470588, -1661.76470588, -1661.76470588,
       -1661.76470588, -1661.76470588,  -833.67083256,  -833.67083256,
        -833.67083256,  -833.67083256,  -833.67083256,  -833.67083256,
        -833.67083256,    -5.57695924,    -5.57695924,    -5.57695924,
          -5.57695924,    -5.57695924,    -5.57695924,    -5.57695924,
         822.51691408,   822.51691408,   822.51691408,   822.51691408,
         822.51691408,   822.51691408,   822.51691408,  1650.6107874 ,
        1650.6107874 ,  1650.6107874 ,  1650.6107874 ,  1650.6107874 ,
        1650.6107874 ,  1650.6107874 ,  2478.70466072,  2478.70466072,
        2478.70466072,  2478.70466072,  2478.70466072,  2478.70466072,
        2478.70466072])
  # objectives(turbine_x, turbine_y, solar, battery)

  nturbs = len(turbine_x)
  turbine_locs = np.zeros((nturbs,2))
  turbine_locs[:,0] = turbine_x
  turbine_locs[:,1] = turbine_y
  place_solar.set_turbine_locs(turbine_locs)
  place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
  place_solar.place_solar()
  solar_poly = place_solar.solar_geometry

  plot_poly(solar_poly,ax=ax9,color="C1")
  plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax9)

  ax9.axis("off")
  ax9.axis("square")


  plt.tight_layout()
  plt.show()