from tkinter import Place
from plotting_functions import plot_turbines, plot_poly
from make_solar import PlaceSolar
from matplotlib import pyplot
import numpy as np


def place_turbines(nturbs,side,rotor_diameter,min_spacing):
    turbine_x = np.array([])
    turbine_y = np.array([])
    for i in range(nturbs):
        placed = False
        while placed == False:
            temp_x = np.random.rand()*side
            temp_y = np.random.rand()*side
            good_point = True
            for j in range(len(turbine_x)):
                dist = np.sqrt((temp_y - turbine_y[j])**2 + (temp_x - turbine_x[j])**2)
                if dist < min_spacing:
                    good_point = False
            if good_point == True:
                turbine_x = np.append(turbine_x,temp_x)
                turbine_y = np.append(turbine_y,temp_y)
                placed = True

    return turbine_x, turbine_y

corner = 2500.0
ngrid = 15 # dunno what this should be, affects how coarse the solar grid is
x = np.linspace(-corner,corner,ngrid)
y = np.linspace(-corner,corner,ngrid)
xlocs = [[i for i in x] for j in y]
ylocs = [[j for i in x] for j in x]
grid_locs = np.zeros((np.shape(xlocs)[0],np.shape(xlocs)[1],2))
grid_locs[:,:,0] = xlocs[:]
grid_locs[:,:,1] = ylocs[:]
grid_locs = np.ndarray.tolist(grid_locs)

min_spacing = 200.0
# place_solar = PlaceSolar(grid_locs, min_spacing,max_arrays=1)
place_solar = PlaceSolar(grid_locs, min_spacing,max_arrays=None)

nturbs = 15
np.random.seed(1)
turbine_x, turbine_y = place_turbines(nturbs,2*corner,100.0,500.0)
print(repr(turbine_x))
print(repr(turbine_y))

# turbine_x = turbine_x - corner
# turbine_y = turbine_y - corner
# turbine_locs = np.zeros((nturbs,2))
# turbine_locs[:,0] = turbine_x
# turbine_locs[:,1] = turbine_y
# place_solar.set_turbine_locs(turbine_locs)
# place_solar.nsolar_cells = ngrid**2

# place_solar.place_solar()