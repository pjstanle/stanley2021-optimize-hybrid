# Copyright 2020 NREL

# Licensed under the Apache License, Version 2.0 (the "License"); you may not
# use this file except in compliance with the License. You may obtain a copy of
# the License at http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
# License for the specific language governing permissions and limitations under
# the License.

# See https://floris.readthedocs.io for documentation

import numpy as np
import matplotlib.pyplot as plt

from shapely.geometry import Polygon
from shapely.geometry import Point

import matplotlib as mpl


def find_lengths(x,y,npoints):
    length = np.zeros(len(x)-1)
    for i in range(npoints):
        length[i] = np.sqrt((x[i+1]-x[i])**2+(y[i+1]-y[i])**2)
    return length


def makeBoundary(start,nturbs,boundary_poly):

    nturbs = int(nturbs)


    xBounds,yBounds = boundary_poly.boundary.coords.xy

    if xBounds[-1] != xBounds[0]:
        xBounds = np.append(xBounds,xBounds[0])
        yBounds = np.append(yBounds,yBounds[0])

    nBounds = len(xBounds)
    lenBound = find_lengths(xBounds,yBounds,len(xBounds)-1)
    circumference = sum(lenBound)
    x = np.zeros(nturbs)
    y = np.zeros(nturbs)

    bound_loc = np.linspace(start,start+circumference-circumference/float(nturbs),nturbs)
    for i in range(nturbs):
        if bound_loc[i] > circumference:
            bound_loc[i] = bound_loc[i]%circumference
        while bound_loc[i] < 0.:
            bound_loc[i] += circumference

    for i in range(nturbs):
        done = False
        for j in range(nBounds):
            if done == False:
                if bound_loc[i] < sum(lenBound[0:j+1]):
                    point_x = xBounds[j] + (xBounds[j+1]-xBounds[j])*(bound_loc[i]-sum(lenBound[0:j]))/lenBound[j]
                    point_y = yBounds[j] + (yBounds[j+1]-yBounds[j])*(bound_loc[i]-sum(lenBound[0:j]))/lenBound[j]
                    done = True
                    x[i] = point_x
                    y[i] = point_y

    return x,y


def discrete_grid(x_spacing,y_spacing,shear,rotation,center_x,center_y,boundary_setback,boundary_poly):
    """
    returns grid turbine layout. Assumes the turbines fill the entire plant area

    Args:
    x_spacing (Float): grid spacing in the unrotated x direction (m)
    y_spacing (Float): grid spacing in the unrotated y direction (m)
    shear (Float): grid shear (rad)
    rotation (Float): grid rotation (rad)
    center_x (Float): the x coordinate of the grid center (m)
    center_y (Float): the y coordinate of the grid center (m)
    boundary_poly (Polygon): a shapely Polygon of the wind plant boundary

    Returns
    return_x (Array(Float)): turbine x locations
    return_y (Array(Float)): turbine y locations
    """

    shrunk_poly = boundary_poly.buffer(-boundary_setback)
    if shrunk_poly.area <= 0:
        return np.array([]), np.array([])
    # create grid
    minx, miny, maxx, maxy = shrunk_poly.bounds
    width = maxx-minx
    height = maxy-miny

    center_point = Point((center_x,center_y))
    poly_to_center = center_point.distance(shrunk_poly.centroid)

    width = np.max([width,poly_to_center])
    height = np.max([height,poly_to_center])
    nrows = int(np.max([width,height])/np.min([x_spacing,y_spacing]))*2 + 1
    ncols = nrows

    xlocs = np.arange(0,ncols)*x_spacing
    ylocs = np.arange(0,nrows)*y_spacing
    row_number = np.arange(0,nrows)

    d = np.array([i for x in xlocs for i in row_number])
    layout_x = np.array([x for x in xlocs for y in ylocs]) + d*y_spacing*np.tan(shear)
    layout_y = np.array([y for x in xlocs for y in ylocs])
    
    # rotate
    rotate_x = np.cos(rotation)*layout_x - np.sin(rotation)*layout_y
    rotate_y = np.sin(rotation)*layout_x + np.cos(rotation)*layout_y

    # move center of grid
    rotate_x = (rotate_x - np.mean(rotate_x)) + center_x
    rotate_y = (rotate_y - np.mean(rotate_y)) + center_y

    # get rid of points outside of boundary polygon
    meets_constraints = np.zeros(len(rotate_x),dtype=bool)
    for i in range(len(rotate_x)):
        pt = Point(rotate_x[i],rotate_y[i])
        if shrunk_poly.contains(pt) or shrunk_poly.touches(pt):
            meets_constraints[i] = True

    # arrange final x,y points
    return_x = rotate_x[meets_constraints]
    return_y = rotate_y[meets_constraints]
   
    return return_x, return_y



if __name__=="__main__":

    start = 230.0
    nTurbs = 50
    boundary_poly = Polygon(([0,0],[10000,-500],[12000,1000],[5000,11000],[-2000,1000]))

    bx,by = makeBoundary(start,nTurbs,boundary_poly)

    nrows = 100
    ncols = 40
    spacing_x = 700
    spacing_y = 500
    shear = 0.1
    rotation = 0.2
    center_x = 4800
    center_y = 6000
    boundary_setback = 1000

    gx,gy = discrete_grid(spacing_x,spacing_y,shear,rotation,center_x,center_y,boundary_setback,boundary_poly)
    
    # side = 10000

    

    # x,y = get_turbine_locs(nrows,ncols,farm_width,farm_height,shear,rotation,center_x,center_y,boundary_mult)

    from plotting_functions import plot_poly
    
    plot_poly(boundary_poly,ax=plt.gca())
    shrunk_poly = boundary_poly.buffer(-boundary_setback)
    plot_poly(shrunk_poly,ax=plt.gca())
    plt.plot(bx,by,"o",color="C0")
    plt.plot(gx,gy,"o",color="C1")
    plt.plot(center_x,center_y,"o",color="black")
    plt.axis("equal")
    plt.show()

    
    