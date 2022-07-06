import numpy as np
from shapely.geometry import Polygon, Point


class SimpleHybrid:
    
    
    def __init__(self):
        # wind
        self.turbine_x = np.array([])
        self.turbine_y = np.array([])

        self.floris_model = None
        self.wind_directions = np.array([])
        self.wind_speeds = np.array([])
        self.wind_frequencies = np.array([])

        # solar
        self.solar_x = np.array([])
        self.solar_y = np.array([])
        self.solar_width = np.array([])
        self.solar_height = np.array([])

        self.solar_capacity_mw = np.array([])
        self.gross_solar_array_aep = np.array([])
        self.solar_array_aep = np.array([])
        self.shadow_multiplier = np.array([])
        

        self.solar_capacity_factor = 0.2 
        self.acres_per_MW = 5.0
        self.shadow_x = 2.0 # rotor diameters of the sigma of shadow loss in x (east west)
        self.shadow_y = 0.5 # rotor diameters of the sigma of shadow loss in y (north south)
        self.shadow_scale = 1.0

        # outputs
        self.solar_aep = 0.0
        self.wind_aep = 0.0
        self.plant_aep = 0.0


    def get_wind_aep(self):
        if len(self.turbine_x) > 0:
            self.floris_model.reinitialize_flow_field(layout_array=(self.turbine_x,self.turbine_y))
            self.wind_aep = self.floris_model.get_farm_AEP(self.wind_directions, self.wind_speeds, self.wind_frequencies, limit_ws=True)
        else:
            self.wind_aep = 0.0


    def get_solar_capacity(self):
        narrays = len(self.solar_x)
        self.solar_capacity_mw = np.zeros(narrays)
        smpa = 4046.86 # square meters per acre
        for i in range(narrays):
            solar_area = self.solar_width[i]*self.solar_height[i]
            self.solar_capacity_mw[i] = solar_area/smpa/self.acres_per_MW


    def get_gross_solar_aep(self):
        narrays = len(self.solar_x)
        self.gross_solar_array_aep = np.zeros(narrays)
        for i in range(narrays):
            self.gross_solar_array_aep[i] = self.solar_capacity_mw[i]*1E6*self.solar_capacity_factor*(365.0*24.0)


    def get_solar_aep(self):

        self.get_solar_capacity()
        self.get_gross_solar_aep()

        narrays = len(self.solar_x)
        nturbs = len(self.turbine_x)
        rotor_diameter = self.floris_model.floris.farm.turbines[0].rotor_diameter

        sigma_x = self.shadow_x*rotor_diameter
        sigma_y = self.shadow_y*rotor_diameter
        
        def gaus2d(sx,sy,dx,dy):
            return np.exp(-(dx**2/(2.0*sx**2)+dy**2/(2.0*sy**2)))

        self.shadow_multiplier = np.ones(narrays)
        for i in range(narrays):
            for j in range(nturbs):
                dx = self.turbine_x[j] - self.solar_x[i]
                dy = self.turbine_y[j] - self.solar_y[i]
                mult = 1.0 - gaus2d(sigma_x,sigma_y,dx,dy)*self.shadow_scale
                self.shadow_multiplier[i] *= mult   

        self.solar_array_aep = self.gross_solar_array_aep*self.shadow_multiplier
        self.solar_aep = np.sum(self.solar_array_aep)
            

    def get_plant_aep(self):

        self.get_wind_aep()
        self.get_solar_aep()
        self.plant_aep = self.wind_aep + self.solar_aep
    

if __name__=="__main__":

    import floris.tools as wfct

    simple = SimpleHybrid()
    simple.floris_model = wfct.floris_interface.FlorisInterface("/Users/astanley/Data/turbines/2_5mw_118d_88h.json")

    simple.turbine_x = [0.0]
    simple.turbine_y = [0.0]
    simple.wind_directions = np.array([270.0])
    simple.wind_speeds = np.array([9.0])
    simple.wind_frequencies = np.array([0.5])

    simple.get_wind_aep()
    print(simple.wind_aep)

    simple.solar_width = np.array([100.0])
    simple.solar_height = np.array([100.0])
    simple.solar_y = np.array([0.0])

    sx = np.linspace(0.0,1000.0,100)
    aep = np.zeros(len(sx))
    for i in range(len(sx)):
        simple.solar_x = np.array([sx[i]])
        simple.get_plant_aep()
        aep[i] = simple.solar_aep
    print("solar_aep: ", simple.solar_aep)
    print("solar_capacity_mw: ", simple.solar_capacity_mw)
    print("gross_solar_array_aep: ", simple.gross_solar_array_aep)
    print("solar_array_aep: ", simple.solar_array_aep)
    print("shadow_multiplier: ", simple.shadow_multiplier)

    import matplotlib.pyplot as plt

    plt.plot(sx,aep)
    # plt.ylim(0.0,max(aep))
    plt.show()

    # def gaus2d(sx,sy,dx,dy):
    #     return np.exp(-(dx**2/(2.0*sx**2)+dy**2/(2.0*sy**2)))

    # from mpl_toolkits.mplot3d import Axes3D
    # from matplotlib import cm
    # 

    # sx = 2.0
    # sy = 0.5
    # X = np.linspace(-6.0,6.0,100)
    # Y = np.linspace(-6.0,6.0,100)
    # X, Y = np.meshgrid(X, Y)

    # Z = gaus2d(sx,sy,X,Y)
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # surf = ax.plot_surface(X, Y, Z,cmap=cm.coolwarm)
    # ax.set_xlabel("X")
    # ax.set_ylabel("Y")
    # plt.show()