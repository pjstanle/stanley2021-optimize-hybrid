import numpy as np
import gfopt.greedy_algorithm

import openmdao.api as om

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
import scipy

import os
from dotenv import load_dotenv

from hybrid.sites import SiteInfo, flatirons_site
# from hybrid.wind_source import WindPlant
from hybrid.hybrid_simulation import HybridSimulation
from hybrid.keys import set_developer_nrel_gov_key
from tools.analysis import create_cost_calculator
from hybrid.grid import Grid

## Main
if __name__ == "__main__":

    # Set API key
    load_dotenv()
    NREL_API_KEY = os.getenv("NREL_API_KEY")
    set_developer_nrel_gov_key(NREL_API_KEY)  # Set this key manually here if you are not setting it using the .env

    # Get resource
    lat = flatirons_site['lat']
    lon = flatirons_site['lon']
    site = SiteInfo(flatirons_site)

    wind_size_mw = 20
    interconnection_size_mw = 5
    solar_size_mw = 25

    technologies = {'solar': solar_size_mw,
                    'wind': wind_size_mw,  # mw system capacity
                    'grid': interconnection_size_mw}
                    
    # Create model
    hybrid_plant = HybridSimulation(technologies, site, interconnect_kw=interconnection_size_mw * 1000)
    hybrid_plant.wake_model = 1

    hybrid_plant.solar.system_capacity_kw = solar_size_mw * 1000
    hybrid_plant.ppa_price = 0.1
    hybrid_plant.setup_cost_calculator(create_cost_calculator(interconnection_size_mw))

    N = np.arange(14)+4
    COE_nom = np.zeros(len(N))
    COE_real = np.zeros(len(N))
    s = 77.0

    # for i in range(len(N)):
    #     n = N[i]
    #     print(n)
        # xarr = np.linspace(0.0,s,n)
        # yarr = np.linspace(0.0,s,n)

    n = 5
    xarr = np.linspace(0.0,s*(n-1),n)
    yarr = np.linspace(0.0,s*(n-1),n)
    turbine_x, turbine_y = np.meshgrid(xarr,yarr)
    layout_x = np.ndarray.flatten(turbine_x)
    layout_y = np.ndarray.flatten(turbine_y)

    hybrid_plant.wind.modify_coordinates(layout_x,layout_y)
    cap = len(layout_x)
    hybrid_plant.setup_cost_calculator(create_cost_calculator(cap))
    hybrid_plant.wind.system_capacity_by_num_turbines(cap*1000)
    # hybrid_plant.grid = Grid(site,cap*1000)

    # hybrid_plant.solar


    hybrid_plant.simulate(25)
    COE_nom = hybrid_plant.lcoe_nom.hybrid
    COE_real = hybrid_plant.lcoe_real.hybrid

    

    print("COE nominal: ", COE_nom)
    print("COE real: ", COE_real)
    print("AEP: ", hybrid_plant.annual_energies.hybrid)
        
    # plt.plot(s/(N-1)/77.0,COE_nom,"o",label="nom")
    # plt.plot(s/(N-1)/77.0,COE_real,"o",label="real")
    # plt.plot(N,COE_nom,"o",label="nom")
    # plt.plot(N,COE_real,"o",label="real")
    # plt.legend()
    # plt.show()



