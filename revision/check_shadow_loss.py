import time
import numpy as np
from shapely.geometry import MultiPolygon
from init_plant import init_hybrid_plant
from make_solar import PlaceSolar
from simple_dispatch import SimpleDispatch
from simple_shadow import SimpleShadow


class HybridPlant():

    def __init__(self, rotor_diameter, hub_height, powercurve_data, solar, shadow, battery_dispatch, site, time_array, solar_kw_per_km2=4000.0):
        
        self.rotor_diameter = rotor_diameter
        self.hub_height = hub_height
        self.powercurve_data = powercurve_data
        self.solar = solar
        self.shadow = shadow
        self.battery_dispatch = battery_dispatch
        self.site = site
        self.time_array = time_array
        self.solar_kw_per_km2 = solar_kw_per_km2
        self.ntime = len(self.time_array)
        
        # user input, updated during optimization
        self.turbine_x = np.array([])
        self.turbine_y = np.array([])
        self.battery_size = 0.0
        self.solar_capacity = 0.0
        


        # internal variables
        dx = self.solar.grid_locs[0][1][0] - grid_locs[0][0][0]
        dy = self.solar.grid_locs[1][0][1] - grid_locs[0][0][1]
        self.solar_cell_area = dx*dy
        self.nturbs = 0
        self.turbine_capacity = np.max(powercurve_data['turbine_powercurve_specification']['turbine_power_output'])

        # mostly internal variables, could be used externally
        self.shadow_losses = np.zeros(self.ntime)
        self.wind_power = np.zeros(self.ntime)
        self.solar_power_no_losses = np.zeros(self.ntime)
        self.solar_power_with_losses = np.zeros(self.ntime)
        self.battery_charged = np.zeros(self.ntime)
        self.battery_discharged = np.zeros(self.ntime)
        self.battery_SOC = np.zeros(self.ntime)

        # outputs
        self.hybrid_power_no_battery = np.zeros(self.ntime)
        self.hybrid_power_with_battery = np.zeros(self.ntime)


    def update_turbine_locs(self, turbine_x, turbine_y):
        self.turbine_x = turbine_x
        self.turbine_y = turbine_y
        self.nturbs = len(turbine_x)


    def evaluate_plant(self):

        print("start solar placement")
        s = time.time()
        turbine_locs = np.zeros((self.nturbs,2))
        turbine_locs[:,0] = self.turbine_x
        turbine_locs[:,1] = self.turbine_y
        self.solar.set_turbine_locs(turbine_locs)
        
        self.solar.nsolar_cells = int(np.floor(self.solar_capacity/(self.solar_kw_per_km2*self.solar_cell_area/1E6)))
        self.solar.place_solar()
        solar_poly = self.solar.solar_geometry

        print("solar placement time: ", time.time()-s)

        print("start power init")
        s = time.time()
        wind_capacity_kw = self.nturbs*self.turbine_capacity
        solar_capacity_kw = solar_poly.area/1E6*self.solar_kw_per_km2
        wind_plant, solar_plant = init_hybrid_plant(self.site, self.rotor_diameter, self.hub_height, self.powercurve_data, wind_capacity_kw, solar_capacity_kw)
        print("power init time: ", time.time()-s)

        wind_plant.modify_coordinates(self.turbine_x,self.turbine_y)
        print("start wind eval")
        s = time.time()
        wind_plant.system_model.execute(0)
        print("wind eval time: ", time.time()-s)

        print("start solar eval")
        s = time.time()
        solar_plant.system_model.execute(0)
        print("solar eval time: ", time.time()-s)

        self.wind_power = wind_plant.system_model.Outputs.gen
        self.solar_power_no_losses = solar_plant.system_model.Outputs.gen

        print("start shadow")
        s = time.time()
        self.shadow.turbine_x = self.turbine_x
        self.shadow.turbine_y = self.turbine_y
        self.shadow.solar_polygons = solar_poly

        for i in range(self.ntime):
            t = self.time_array[i]%24.0
            d = self.time_array[i]/24.0
            self.shadow_losses[i] = self.shadow.calculate_losses(t,d)
            self.solar_power_with_losses[i] = self.solar_power_no_losses[i]*(1-self.shadow_losses[i])
        print("shadow eval time: ", time.time()-s)

        print("start battery dispatch")
        s = time.time()
        self.hybrid_power_no_battery = self.wind_power + self.solar_power_with_losses

        self.battery_dispatch.plant_power = self.hybrid_power_no_battery
        self.battery_dispatch.battery_size = self.battery_size

        self.battery_charged, self.battery_discharged, self.battery_SOC = self.battery_dispatch.run()
        self.hybrid_power_with_battery = self.hybrid_power_no_battery - self.battery_charged + self.battery_discharged
        print("battery dispatch eval time: ", time.time()-s)

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

    start_time = 0
    ntime = 8760
    data = site.wind_resource.data["data"] 

    new_data = np.zeros((ntime,4))
    for i in range(ntime):
        new_data[i,0] = data[start_time+i][0]
        new_data[i,1] = data[start_time+i][1]
        new_data[i,2] = data[start_time+i][2]
        new_data[i,3] = data[start_time+i][3]
    
    site.wind_resource.data["data"] = np.ndarray.tolist(new_data)

    corner = 5000
    # make solar
    ngrid = 50
    x = np.linspace(0,corner,ngrid)
    y = np.linspace(0,corner,ngrid)
    xlocs = [[i for i in x] for j in y]
    ylocs = [[j for i in x] for j in x]
    grid_locs = np.zeros((np.shape(xlocs)[0],np.shape(xlocs)[1],2))
    grid_locs[:,:,0] = xlocs[:]
    grid_locs[:,:,1] = ylocs[:]
    grid_locs = np.ndarray.tolist(grid_locs)
    min_spacing = rotor_diameter*1.0
    solar = PlaceSolar(grid_locs,min_spacing)

    # shadow
    tower_width = 5.0
    shadow = SimpleShadow(hub_height, rotor_diameter, tower_width, lat, lon)

    # battery
    min_power = 1000.0
    battery_dispatch = SimpleDispatch(ntime, min_power)

    time_array = np.arange(ntime) + start_time
    hybrid_plant = HybridPlant(rotor_diameter, hub_height, powercurve_data, solar, shadow, battery_dispatch, site, time_array)

    turbine_x = np.array([0,corner,corner])
    turbine_y = np.array([0,0,corner])
    battery_size = 10000
    solar_capacity = 15000


    hybrid_plant.update_turbine_locs(turbine_x, turbine_y)
    hybrid_plant.battery_size = battery_size
    hybrid_plant.solar_capacity = solar_capacity

    hybrid_plant.evaluate_plant()


    shadow_losses = hybrid_plant.shadow_losses
    wind_power = np.array(hybrid_plant.wind_power)
    solar_power_no_losses = np.array(hybrid_plant.solar_power_no_losses)
    solar_power_with_losses = hybrid_plant.solar_power_with_losses
    battery_charged = hybrid_plant.battery_charged
    battery_discharged = hybrid_plant.battery_discharged
    battery_SOC = hybrid_plant.battery_SOC
    hybrid_power_no_battery = hybrid_plant.hybrid_power_no_battery
    hybrid_power_with_battery = hybrid_plant.hybrid_power_with_battery

    plt.figure(1)
    plt.plot(shadow_losses,label="shadow losses")
    plt.plot(solar_power_no_losses/max(solar_power_no_losses), label="no losses")
    plt.plot(solar_power_with_losses/max(solar_power_no_losses), label="with losses")
    plt.legend()
    plt.title("shadow_losses")

    plt.figure(2)
    plt.plot(solar_power_with_losses,label="solar power")
    plt.plot(wind_power, label="wind power")
    plt.plot(hybrid_power_no_battery, label="hybrid power")
    plt.legend()
    plt.title("plant power")

    plt.figure(3)
    plt.plot(battery_charged, label="charge")
    plt.plot(battery_discharged, label="discharge")
    plt.plot(battery_SOC, label="SOC")
    plt.legend()
    plt.title("battery")

    plt.figure(4)
    plt.plot(hybrid_power_no_battery, label="without battery")
    plt.plot(battery_charged, label="charge")
    plt.plot(battery_discharged, label="discharge")
    plt.plot(hybrid_power_with_battery, label="with battery")
    plt.legend()
    plt.title("total power")

    plt.show()