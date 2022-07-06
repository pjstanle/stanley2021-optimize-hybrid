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


import os

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

    
def makeGrid(nrows,ncols,spacing_x,spacing_y,shear,rotation,center_x,center_y,boundary_setback,boundary_poly):

    # create grid
    nrows = int(nrows)
    ncols = int(ncols)
    xlocs = np.arange(ncols)*spacing_x
    ylocs = np.arange(nrows)*spacing_y

    nturbs = nrows*ncols
    layout_x = np.zeros(nturbs)
    layout_y = np.zeros(nturbs)
    turb = 0
    for i in range(nrows):
        for j in range(ncols):
            layout_x[turb] = xlocs[j] + float(i)*spacing_y*np.tan(shear)
            layout_y[turb] = ylocs[i]
            turb += 1
    
    # rotate
    rotate_x = np.cos(rotation)*layout_x - np.sin(rotation)*layout_y
    rotate_y = np.sin(rotation)*layout_x + np.cos(rotation)*layout_y

    # move center of grid
    centered_x = (rotate_x - np.mean(rotate_x)) + center_x
    centered_y = (rotate_y - np.mean(rotate_y)) + center_y

    # get rid of points outside of boundary
    turbs = np.zeros((len(centered_x),2))
    turbs[:,0] = centered_x[:]
    turbs[:,1] = centered_y[:]

    grid_poly = boundary_poly.buffer(-boundary_setback)
    bgx,bgy = grid_poly.boundary.coords.xy
    edges = np.zeros((len(bgx),2))
    edges[:,0] = bgx[:]
    edges[:,1] = bgy[:]
    boundary = mpl.path.Path(edges)
    in_bounds = boundary.contains_points(turbs)

    return_x = np.zeros(sum(in_bounds))
    return_y = np.zeros(sum(in_bounds))
    ind = 0
    for i in range(len(rotate_x)):
        if in_bounds[i]:
            return_x[ind] = centered_x[i]
            return_y[ind] = centered_y[i]
            ind += 1

    return return_x, return_y



if __name__=="__main__":

    start = 230.0
    nTurbs = 10
    boundary_poly = Polygon(([0,0],[1000,-1000],[1200,1000],[500,1100],[-200,1000])).buffer(2000)

    bx,by = makeBoundary(start,nTurbs,boundary_poly)

    nrows = 30
    ncols = 40
    spacing_x = 1500
    spacing_y = 900
    shear = 0.1
    rotation = 0.2
    center_x = 5000
    center_y = 6000
    boundary_setback = 700

    gx,gy = makeGrid(nrows,ncols,spacing_x,spacing_y,shear,rotation,center_x,center_y,boundary_setback,boundary_poly)
    
    # side = 10000

    rotor_diameter = 200.0

    # x,y = get_turbine_locs(nrows,ncols,farm_width,farm_height,shear,rotation,center_x,center_y,boundary_mult)

    ox,oy = boundary_poly.exterior.coords.xy

    plt.figure(figsize=(4,4))
    plt.plot(ox,oy,"--k",linewidth=0.5)
    turbine_x = np.append(gx,bx)
    turbine_y = np.append(gy,by)
    for i in range(len(turbine_x)):
        turb = plt.Circle((turbine_x[i],turbine_y[i]),rotor_diameter/2.0,facecolor="C0")
        plt.gca().add_patch(turb)
    plt.axis("equal")
    plt.axis("off")

    plt.subplots_adjust(top=0.99,bottom=0.01,left=0.01,right=0.99)

    plt.savefig("figures/bg_turbines.pdf",transparent=True)
    plt.show()

    
    