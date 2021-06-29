import numpy as np
from shapely.geometry import Polygon, Point, LineString, MultiPolygon
from gradient_free import GeneticAlgorithm
from model_hybrid import SimpleHybrid
import floris.tools as wfct
import matplotlib.pyplot as plt
import time
from shapely.ops import unary_union
from boundary_grid import makeBoundary, makeGrid
from simple_flicker import SimpleFlicker


def calculate_flicker(turbine_x,turbine_y,solar_poly):

        global rotor_diameter
        global plant

        def _gaus2d(sx,sy,dx,dy):
            return np.exp(-(dx**2/(2.0*sx**2)+dy**2/(2.0*sy**2)))

        nturbs = len(turbine_x)

        sigma_x = plant.shadow_x*rotor_diameter
        sigma_y = plant.shadow_y*rotor_diameter
        solar_loc = solar_poly.centroid
        shadow_multiplier = 1.0
        for j in range(nturbs):
            dx = turbine_x[j] - solar_loc.x
            dy = turbine_y[j] - solar_loc.y
            mult = 1.0 - _gaus2d(sigma_x,sigma_y,dx,dy)
            shadow_multiplier *= mult 
        
        return shadow_multiplier


def get_solar_polygon(solar_x,solar_y,solar_width,solar_height):

    solar_pts = [(solar_x-solar_width/2.0,solar_y-solar_height/2.0),
                (solar_x+solar_width/2.0,solar_y-solar_height/2.0),
                (solar_x+solar_width/2.0,solar_y+solar_height/2.0),
                (solar_x-solar_width/2.0,solar_y+solar_height/2.0),
                (solar_x-solar_width/2.0,solar_y-solar_height/2.0)]

    
    solar_poly = Polygon(solar_pts)
    return solar_poly


def get_wind_solar_distance(turbine_x,turbine_y,solar_poly):
              
    global rotor_diameter
    
    x,y = solar_poly.boundary.coords.xy
    solar_pts = np.zeros((len(x),2))
    solar_pts[:,0] = x[:]
    solar_pts[:,1] = y[:]

    solar_line = LineString(solar_pts)
    turbine_point = Point(turbine_x,turbine_y)
    spacing = solar_line.distance(turbine_point)
    if solar_poly.contains(turbine_point):
        spacing *= -1.0

    return spacing-rotor_diameter


def plot_hybrid(turbine_x,turbine_y,solar_polys):

    global boundary_poly
    global rotor_diameter

    plt.cla()

    # plot boundaries
    bx,by = boundary_poly.boundary.coords.xy
    bx = np.append(bx,bx[0])
    by = np.append(by,by[0])
    plt.plot(bx,by,"--k",linewidth=0.5)

    # plot solar
    try:
        for j in range(len(solar_polys)):
            x,y = solar_polys[j].boundary.coords.xy
            plt.fill(x,y,color="C1")
    except:
        x,y = solar_polys.boundary.coords.xy
        plt.fill(x,y,color="C1")

    # plot turbines
    for i in range(len(turbine_x)):
        turb = plt.Circle((turbine_x[i],turbine_y[i]),rotor_diameter/2.0,facecolor="C0")
        plt.gca().add_patch(turb)

    
    plt.axis("equal")
    plt.pause(0.0001)


if __name__=="__main__":

    global rotor_diameter
    global plant
    global boundary_poly


    best_answer = 1E20

    # operation = "always_stay_charged"
    operation = "store_curtailment"
    solar_array_penalty = 100.0
    function_calls = 0

    battery_cost = 1.0
    under_penalty = 0.5
    solar_cost = 1.0
    min_energy = 40.0
    interconnect = 50.0

    boundary_poly = Polygon(([0,0],[1000,-100],[1200,1000],[500,1100],[-200,1000])).buffer(1000)

    np.random.seed(seed=1)

    ws_array = np.array([0.95,1.00,1.00,0.95,0.85,0.7,0.5,0.3,0.15,0.05,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.15,0.3,0.5,0.7,0.85])*6.0+6.0
    sr_array = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.15,0.3,0.5,0.7,0.85,0.95,1.00,1.00,0.95,0.85,0.7,0.5,0.3,0.15,0.05,0.0,0.0])
    ws_array = np.append(ws_array,ws_array)
    sr_array = np.append(sr_array,sr_array)

    ws_array = ws_array - np.random.rand(len(ws_array))
    sr_array = sr_array - np.random.rand(len(sr_array))*0.1
    for i in range(len(sr_array)):
        if sr_array[i] < 0.0:
            sr_array[i] = 0.0

    wd_array = np.ones(len(ws_array))*270.0 + np.random.rand(len(ws_array))*30.0
    
    plant = SimpleHybrid()
    plant.floris_model = wfct.floris_interface.FlorisInterface("/Users/astanley/Data/turbines/2_5mw_118d_88h_1jensen.json")
    rotor_diameter = plant.floris_model.floris.farm.turbines[0].rotor_diameter
    plant.acres_per_MW = 5.0
    plant.shadow_x = 4.0 # rotor diameters of the sigma of shadow loss in x (east west)
    plant.shadow_y = 0.5 # rotor diameters of the sigma of shadow loss in y (north south)

    D = 118.0
    wind_wind_spacing = D*4.0
    wind_solar_spacing = D*0.5

    # turbine variables
    boundary_start = 0.0
    boundary_turbs = 5
    grid_rows = 20
    grid_cols = 20
    grid_spacing_x = 1000
    grid_spacing_y = 1000
    grid_shear = 0
    grid_rotation = 0
    grid_center_x = 500
    grid_center_y = 500
    grid_boundary_setback = 200

    # solar variables
    solar_x = [-500,0,500]
    solar_y = [-500,0,500]
    solar_width = [400,400,400]
    solar_height = [400,400,400]

    # assign turbine variables
    boundary_turbine_x, boundary_turbine_y = makeBoundary(boundary_start,boundary_turbs,boundary_poly)
    grid_turbine_x, grid_turbine_y = makeGrid(grid_rows,grid_cols,grid_spacing_x,grid_spacing_y,grid_shear,
                                                grid_rotation,grid_center_x,grid_center_y,grid_boundary_setback,boundary_poly)
    initial_turbine_x = np.append(boundary_turbine_x,grid_turbine_x)
    initial_turbine_y = np.append(boundary_turbine_y,grid_turbine_y)

    poly0 = get_solar_polygon(solar_x[0],solar_y[0],solar_width[0],solar_height[0]).intersection(boundary_poly)
    poly1 = get_solar_polygon(solar_x[1],solar_y[1],solar_width[1],solar_height[1]).intersection(boundary_poly)
    poly2 = get_solar_polygon(solar_x[2],solar_y[2],solar_width[2],solar_height[2]).intersection(boundary_poly)

    nsolar = 0
    if poly0.area > 0.0:
        nsolar += 1
    if poly1.area > 0.0:
        nsolar += 1
    if poly2.area > 0.0:
        nsolar += 1
    
    solar_polys = [0]*nsolar
    counter = 0
    if poly0.area > 0.0:
        solar_polys[counter] = poly0
        counter += 1
    if poly1.area > 0.0:
        solar_polys[counter] = poly1
        counter += 1
    if poly2.area > 0.0:
        solar_polys[counter] = poly2
        counter += 1
    
    polygons = [Point(i, 0).buffer(0.7) for i in range(5)]
    solar_union = unary_union(polygons)
    solar_union = unary_union(solar_polys)
    if type(solar_union) == Polygon:
        solar_union = MultiPolygon([solar_union])


    # delete turbines close to solar
    keep_turbines = np.ones(len(initial_turbine_x),dtype=int)
    for i in range(len(initial_turbine_x)):
        for j in range(len(solar_union)):
            distance = get_wind_solar_distance(initial_turbine_x[i],initial_turbine_y[i],solar_union[j])
            if distance < wind_solar_spacing:
                keep_turbines[i] = 0
    
    turbine_x = np.zeros(np.sum(keep_turbines))
    turbine_y = np.zeros(np.sum(keep_turbines))
    counter = 0
    for i in range(len(initial_turbine_x)):
        if keep_turbines[i] == 1:
            turbine_x[counter] = initial_turbine_x[i]
            turbine_y[counter] = initial_turbine_y[i]
            counter += 1

    shadow_scale = 1.0
    loss0 = 1-calculate_flicker(turbine_x,turbine_y,poly0)**shadow_scale
    loss1 = 1-calculate_flicker(turbine_x,turbine_y,poly1)**shadow_scale
    loss2 = 1-calculate_flicker(turbine_x,turbine_y,poly2)**shadow_scale

    print("Gaus")
    print(loss0)
    print(loss1)
    print(loss2)

    solar_x,solar_y = solar_union[0].boundary.coords.xy
    T = 15.0
    flicker = SimpleFlicker()
    flicker.turbine_x = turbine_x
    flicker.turbine_y = turbine_y
    flicker.solar_x = solar_x
    flicker.solar_y = solar_y
    flicker.calculate_losses(T)

    # flicker.lurbine_locs = turbine_locs
    # class SimpleFlicker():

    # def __init__(self, solar_verts, T, turbine_locs):

    #     self.turbine_locs = [[0, 0]]
    #     self.solar_verts = solar_verts
    #     self.turbine_locs = turbine_locs

    plot_hybrid(turbine_x,turbine_y,solar_union)
    plt.show()