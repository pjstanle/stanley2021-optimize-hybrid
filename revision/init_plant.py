from hybrid.wind_source import WindPlant
from hybrid.solar_source import SolarPlant
from hybrid.log import hybrid_logger as logger

import numpy as np


def init_wind_plant(site,tower_height,rotor_diameter,powercurve_data,system_capacity_kw):

    # Set up HOPP models
    wind_plant = WindPlant(site, system_capacity_kw, rotor_diameter, tower_height, powercurve_data)

    wind_plant.system_model.Turbine.wind_turbine_hub_ht = tower_height
    wind_plant.system_model.Turbine.wind_turbine_rotor_diameter = rotor_diameter
    wind_plant.system_model.Turbine.wind_turbine_powercurve_windspeeds = \
        powercurve_data['turbine_powercurve_specification']['wind_speed_ms']
    wind_plant.system_model.Turbine.wind_turbine_powercurve_powerout = \
        powercurve_data['turbine_powercurve_specification']['turbine_power_output']

    wind_plant.wake_model = 1


    # set time carying losses to zero
    wind_plant.system_model.Losses.avail_bop_loss = 0.0
    wind_plant.system_model.Losses.avail_grid_loss = 0.0
    wind_plant.system_model.Losses.avail_turb_loss = 0.0
    wind_plant.system_model.Losses.elec_eff_loss = 0.0
    wind_plant.system_model.Losses.elec_parasitic_loss = 0.0
    wind_plant.system_model.Losses.env_degrad_loss = 0.0
    wind_plant.system_model.Losses.env_env_loss = 0.0
    wind_plant.system_model.Losses.env_exposure_loss = 0.0
    wind_plant.system_model.Losses.ops_env_loss = 0.0
    wind_plant.system_model.Losses.ops_grid_loss = 0.0
    wind_plant.system_model.Losses.ops_load_loss = 0.0
    wind_plant.system_model.Losses.ops_strategies_loss = 0.0
    wind_plant.system_model.Losses.turb_generic_loss = 0.0
    wind_plant.system_model.Losses.turb_hysteresis_loss = 0.0
    wind_plant.system_model.Losses.turb_perf_loss = 0.0
    wind_plant.system_model.Losses.turb_specific_loss = 0.0
    wind_plant.system_model.Losses.wake_ext_loss = 0.0
    wind_plant.system_model.Losses.wake_future_loss = 0.0
    wind_plant.system_model.Losses.wake_int_loss = 0.0


    # print("assign: ", wind_plant.system_model.Losses.assign)
    # print("avail_bop_loss", wind_plant.system_model.Losses.avail_bop_loss)
    # print("avail_grid_loss: ", wind_plant.system_model.Losses.avail_grid_loss)
    # print("avail_turb_loss: ", wind_plant.system_model.Losses.avail_turb_loss)
    # print("elec_eff_loss", wind_plant.system_model.Losses.elec_eff_loss)
    # print("elec_parasitic_loss: ", wind_plant.system_model.Losses.elec_parasitic_loss)
    # # print("en_icing_cutoff: ", wind_plant.system_model.Losses.en_icing_cutoff)
    # # print("en_low_temp_cutoff: ", wind_plant.system_model.Losses.en_low_temp_cutoff)
    # print("env_degrad_loss: ", wind_plant.system_model.Losses.env_degrad_loss)
    # print("env_env_loss: ", wind_plant.system_model.Losses.env_env_loss)
    # print("env_exposure_loss: ", wind_plant.system_model.Losses.env_exposure_loss)
    # print("export: ", wind_plant.system_model.Losses.export)
    # # print("icing_cutoff_rh: ", wind_plant.system_model.Losses.icing_cutoff_rh)
    # # print("icing_cutoff_temp: ", wind_plant.system_model.Losses.icing_cutoff_temp)
    # # print("low_temp_cutoff: ", wind_plant.system_model.Losses.low_temp_cutoff)
    # print("ops_env_loss: ", wind_plant.system_model.Losses.ops_env_loss)
    # print("ops_grid_loss: ", wind_plant.system_model.Losses.ops_grid_loss)
    # print("ops_load_loss: ", wind_plant.system_model.Losses.ops_load_loss)
    # print("ops_strategies_loss: ", wind_plant.system_model.Losses.ops_strategies_loss)
    # print("replace: ", wind_plant.system_model.Losses.replace)
    # print("turb_generic_loss: ", wind_plant.system_model.Losses.turb_generic_loss)
    # print("turb_hysteresis_loss: ", wind_plant.system_model.Losses.turb_hysteresis_loss)
    # print("turb_perf_loss: ", wind_plant.system_model.Losses.turb_perf_loss)
    # print("turb_specific_loss: ", wind_plant.system_model.Losses.turb_specific_loss)
    # print("wake_ext_loss: ", wind_plant.system_model.Losses.wake_ext_loss)
    # print("wake_future_loss: ", wind_plant.system_model.Losses.wake_future_loss)
    # print("wake_int_loss: ", wind_plant.system_model.Losses.wake_int_loss)
    # print(dir(wind_plant.system_model.Losses))

    logger.disabled = True

    return wind_plant


def init_solar_plant(site,system_capacity_kw):

    solar_plant = SolarPlant(site, system_capacity_kw)

    return solar_plant


def init_hybrid_plant(site, rotor_diameter, tower_height, powercurve_data, wind_capacity_kw, solar_capacity_kw):
    wind_plant = init_wind_plant(site,tower_height,rotor_diameter, powercurve_data, wind_capacity_kw)
    solar_plant = init_solar_plant(site, solar_capacity_kw)

    return wind_plant, solar_plant


if __name__=="__main__":
    import os
    from dotenv import load_dotenv
    from hybrid.sites import SiteInfo
    from hybrid.sites import flatirons_site as sample_site
    from hybrid.keys import set_developer_nrel_gov_key
    import json
    import matplotlib.pyplot as plt

    # Set API key
    load_dotenv()
    NREL_API_KEY = os.getenv("NREL_API_KEY")
    set_developer_nrel_gov_key(NREL_API_KEY)  # Set this key manually here if you are not setting it using the .env
    
    powercurve_filename = 'turbine_data/high_7r_200d_135h.txt'
    rotor_diameter = 200.0
    hub_height = 135.0
    turbine_rating = 7.0

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
    # site = SiteInfo(sample_site)

    data = site.wind_resource.data["data"]
    nhours = 8760
    for k in range(nhours):
        data[k][2] = data[k][2]
    site.wind_resource.data["data"] = data   

    

    print("WIND")
    N = 10
    np.random.seed(0)
    turbine_x = np.random.rand(N)*10000
    turbine_y = np.random.rand(N)*10000
    wind_plant = init_wind_plant(site,hub_height,rotor_diameter)
    wind_plant.modify_coordinates(turbine_x,turbine_y)
    wind_plant.simulate(1)
    print(wind_plant.annual_energy_kw())
    print(wind_plant.annual_energy_kw()/N)
    print(wind_plant.system_capacity_kw)

    N = 30
    turbine_x = np.random.rand(N)*10000
    turbine_y = np.random.rand(N)*10000
    wind_plant = init_wind_plant(site,hub_height,rotor_diameter)
    wind_plant.modify_coordinates(turbine_x,turbine_y)
    wind_plant.simulate(1)
    print(wind_plant.annual_energy_kw())
    print(wind_plant.annual_energy_kw()/N)
    print(wind_plant.system_capacity_kw)


    print("SOLAR")

    solar_plant = init_solar_plant(site)
    solar_plant.system_capacity_kw = 10000    
    solar_plant.simulate(1)
    print(solar_plant.annual_energy_kw())

    solar_plant = init_solar_plant(site)
    solar_plant.system_capacity_kw = 20000 
    solar_plant.simulate(1)
    print(solar_plant.annual_energy_kw())

    # print(wind_plant.system_model.Outputs)
    gen = wind_plant.system_model.Outputs.gen
    print(gen)
    plt.plot(gen)
    plt.show()    