import os
from dotenv import load_dotenv

from hybrid.sites import SiteInfo, flatirons_site
from hybrid.wind_source import WindPlant
from hybrid.keys import set_developer_nrel_gov_key

import numpy as np


# Set API key
load_dotenv()
NREL_API_KEY = os.getenv("NREL_API_KEY")

print("key: ", NREL_API_KEY)

set_developer_nrel_gov_key(NREL_API_KEY)  # Set this key manually here if you are not setting it using the .env

# Set wind, solar, and interconnection capacities (in MW)
wind_size_mw = 20

# Get resource
lat = flatirons_site['lat']
lon = flatirons_site['lon']
site = SiteInfo(flatirons_site)

# Create model
wind_plant = WindPlant(site, wind_size_mw)
wind_plant.wake_model = 1
print(wind_plant.wake_model)

nx = 10
sx = 600.0
ny = 10
sy = 600.0
xarr = np.linspace(0.0,(nx-1)*sx,nx)
yarr = np.linspace(0.0,(ny-1)*sy,ny)
turbine_x, turbine_y = np.meshgrid(xarr,yarr)
turbine_x = np.ndarray.flatten(turbine_x)
turbine_y = np.ndarray.flatten(turbine_y)

wind_plant.modify_coordinates(turbine_x,turbine_y)
# wind_plant.system_capacity_by_num_turbines(19.5 * 1000)

print("num turbines: ", wind_plant.num_turbines)
print("system capacity: ", wind_plant.system_capacity_kw)
print("turb rating: ", wind_plant.turb_rating)

wind_plant.simulate(1)
print("AEP: ", wind_plant.annual_energy_kw())
print("rotor diameter: ", wind_plant.rotor_diameter)