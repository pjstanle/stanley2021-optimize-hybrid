import numpy as np
import gfopt.greedy_algorithm

import os
from dotenv import load_dotenv

from hybrid.sites import SiteInfo, flatirons_site
from hybrid.solar_source import SolarPlant
from hybrid.keys import set_developer_nrel_gov_key

import PySAM.Pvwattsv7 as pvwatts
from hybrid.flicker.flicker_mismatch import module_width, module_height, FlickerMismatch

import numpy as np
import matplotlib.pyplot as plt

## Main
if __name__ == "__main__":

    # Set API key
    load_dotenv(    )
    NREL_API_KEY = os.getenv("NREL_API_KEY")
    set_developer_nrel_gov_key(NREL_API_KEY)  # Set this key manually here if you are not setting it using the .env

    # Set wind, solar, and interconnection capacities (in MW)
    solar_size_mw = 20

    # Get resource
    lat = flatirons_site['lat']
    lon = flatirons_site['lon']
    site = SiteInfo(flatirons_site)

    # Create model

    solar_size_array = np.linspace(0.0001,1000.0,1000)
    aep = np.zeros(len(solar_size_array))
    for i in range(len(solar_size_array)):
        print(i)
        solar_plant = SolarPlant(site, solar_size_array[i])
        solar_plant.simulate(1)
        aep[i] = solar_plant.annual_energy_kw()

    print(repr(aep))
    plt.plot(solar_size_array,aep)
    plt.show()

    
