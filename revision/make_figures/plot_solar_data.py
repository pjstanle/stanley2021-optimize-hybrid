
import numpy as np
import matplotlib.pyplot as plt
import os
from dotenv import load_dotenv
from hybrid.sites import SiteInfo
from hybrid.sites import flatirons_site as sample_site
from hybrid.keys import set_developer_nrel_gov_key
import scipy.interpolate


if __name__=="__main__":



        
    load_dotenv()
    NREL_API_KEY = os.getenv("NREL_API_KEY")
    set_developer_nrel_gov_key(NREL_API_KEY)  # Set this key manually here if you are not setting it using the .env
    lat = 40.966
    lon = -84.598
    year = 2012
    sample_site['year'] = year
    sample_site['lat'] = lat
    sample_site['lon'] = lon
    hub_height = 88.0
    site = SiteInfo(sample_site, hub_height=hub_height)



    plt.figure(figsize=(5,2.25))

    ax = plt.subplot(121)
    ax.set_title("cost curve",fontsize=8)
    # capex_cost_real = np.array([2*1382.0,1382.0,1124.0,966.0,887.0,849.0,792.0,765.0]) # $/kW realistic
    # capex_size = np.array([1.0,20.0,50.0,100.0,150.0,200.0,400.0,1000.0]) # MW
    # cost = capex_size*capex_cost_real*1000.0
    # f = scipy.interpolate.interp1d(capex_size,cost,kind='cubic')

    solar_cost = np.array([5*1786.0,1786.0,1622.0,1528.0,1494.0,1470.0,1421.0,1408.0])*1555.0/1494.0 *1075/1555 # $/kW/year realistic
    solar_capacity = np.array([0.0,20.0,50.0,100.0,150.0,200.0,400.0,1000.0]) # MW
    f = scipy.interpolate.interp1d(solar_capacity, solar_cost*solar_capacity, kind='cubic')
    
    x = np.linspace(1E-6,300.0,100)
    ax.plot(x,f(x)/x,color="C2")

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis='both', labelsize=8)
    ax.set_xlabel("solar capacity (MW)",fontsize=8)

    # ax.set_ylim(1550.0,max(f(x)/x))
    # ax.set_xlim((0,600))

    # ax.set_yticks((0.8,1.0,1.2,1.4,1.6,1.8))
    ax.set_ylabel("capex ($M/kW)",fontsize=8)
    # ax.set_yticklabels(("0.8","1.0","1.2","1.4","1.6","1.8"))

    ax2 = plt.subplot(122)
    # data = np.array(site.wind_resource.data["data"])
    data = site.solar_resource.data
    print(data.keys())

    plt.figure(2)
    plt.plot(data["dn"])
    plt.figure(3)
    plt.plot(data["df"])
    plt.figure(4)
    plt.plot(data["gh"])
    # speed = data[:,2]
    # direction = data[:,3]
    # hours = np.arange(8760)

    # start = 0
    # nhours = 72
    
    # plt.plot(hours[start:start+nhours],direction[start:start+nhours])
    # plt.xlabel("time (hr)",fontsize=8)
    # plt.ylabel("wind direction (deg)",fontsize=8)
    # plt.gca().set_xticks((start+0,start+12,start+24,start+36,start+48,start+60,start+72))
    # plt.gca().set_xticklabels(("0","12","24","36","48","60","72"))
    # plt.xticks(fontsize=8)
    # plt.yticks(fontsize=8)

    # plt.subplot(122)
    # plt.plot(hours[start:start+nhours],speed[start:start+nhours],color="C1")
    # plt.xlabel("time (hr)",fontsize=8)
    # plt.ylabel("wind speed (m/s)",fontsize=8)
    # plt.gca().set_xticks((start+0,start+12,start+24,start+36,start+48,start+60,start+72))
    # plt.gca().set_xticklabels(("0","12","24","36","48","60","72"))
    # plt.xticks(fontsize=8)
    # plt.yticks(fontsize=8)

    plt.tight_layout()
    # # plt.savefig("hourly_wind_resource.pdf", transparent=True)
    plt.show()