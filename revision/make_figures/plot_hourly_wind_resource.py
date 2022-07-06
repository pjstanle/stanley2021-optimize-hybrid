
import numpy as np
import matplotlib.pyplot as plt
import os
from dotenv import load_dotenv
from hybrid.sites import SiteInfo
from hybrid.sites import flatirons_site as sample_site
from hybrid.keys import set_developer_nrel_gov_key


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
    data = np.array(site.wind_resource.data["data"])

    speed = data[:,2]
    direction = data[:,3]
    hours = np.arange(8760)

    start = 0
    nhours = 72
    plt.figure(figsize=(5,2.25))
    plt.subplot(121)
    plt.plot(hours[start:start+nhours],direction[start:start+nhours])
    plt.xlabel("time (hr)",fontsize=8)
    plt.ylabel("wind direction (deg)",fontsize=8)
    plt.gca().set_xticks((start+0,start+12,start+24,start+36,start+48,start+60,start+72))
    plt.gca().set_xticklabels(("0","12","24","36","48","60","72"))
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)

    plt.subplot(122)
    plt.plot(hours[start:start+nhours],speed[start:start+nhours],color="C1")
    plt.xlabel("time (hr)",fontsize=8)
    plt.ylabel("wind speed (m/s)",fontsize=8)
    plt.gca().set_xticks((start+0,start+12,start+24,start+36,start+48,start+60,start+72))
    plt.gca().set_xticklabels(("0","12","24","36","48","60","72"))
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)

    plt.tight_layout()
    # plt.savefig("hourly_wind_resource.pdf", transparent=True)
    plt.show()