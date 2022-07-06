from matplotlib.pyplot import grid
import numpy as np
import shapely.geometry as geom
import time

from shapely.geometry.polygon import Polygon, Point, LineString


class PlaceSolar():
    
    def __init__(self, grid_locs, min_spacing, max_arrays=None):

        self.grid_locs = grid_locs

        self.turbine_locs = []
        self.nsolar_cells = 0
        if max_arrays == None:
            self.max_arrays = 1E16
        else:
            self.max_arrays = max_arrays
        
        self.min_spacing = min_spacing

        self.nturbs = len(self.turbine_locs)
        self.dx = grid_locs[0][1][0] - grid_locs[0][0][0]
        self.dy = grid_locs[1][0][1] - grid_locs[0][0][1]
        self.nrows, self.ncols, _ = np.shape(grid_locs)
        self.distance = np.zeros((self.nrows,self.ncols))
        self.solar_indices = []
        self.solar_geometry = geom.Polygon()
        
        self.solar_locs = []


    def reset(self):
        self.distance = np.zeros((self.nrows,self.ncols))
        self.solar_indices = []
        self.solar_geometry = geom.Polygon()


    def set_turbine_locs(self, turbine_locs):
        self.turbine_locs = turbine_locs
        self.nturbs = len(self.turbine_locs)


    def find_distances(self):

        turbine_points = np.zeros(self.nturbs, dtype=Point)
        for i in range(self.nturbs):
            turbine_points[i] = Point((self.turbine_locs[i][0],self.turbine_locs[i][1]))

        for i in range(self.nrows):
            for j in range(self.ncols):
                solar_cell = Polygon(([self.grid_locs[i][j][0]-self.dx/2.0, self.grid_locs[i][j][1]-self.dy/2.0],
                                      [self.grid_locs[i][j][0]-self.dx/2.0, self.grid_locs[i][j][1]+self.dy/2.0],
                                      [self.grid_locs[i][j][0]+self.dx/2.0, self.grid_locs[i][j][1]+self.dy/2.0],
                                      [self.grid_locs[i][j][0]+self.dx/2.0, self.grid_locs[i][j][1]-self.dy/2.0]))

                d = np.zeros(self.nturbs)
                for k in range(self.nturbs):
                    d[k] = turbine_points[k].distance(solar_cell)
                self.distance[i,j] = np.min(d)

    def find_unattached_cell(self):
        index = np.unravel_index(self.distance.argmax(), self.distance.shape)
        cell_index = [index[0],index[1]]
        cell_location = self.grid_locs[cell_index[0]][cell_index[1]]
        self.distance[cell_index[0]][cell_index[1]] = -1E6
        
        return cell_location, cell_index
        
    def find_attached_cell(self):
        attached_indices = []
        for i in range(len(self.solar_indices)):
            xind, yind = self.solar_indices[i]
            if yind + 1 <= self.nrows - 1 and [xind, yind + 1] not in attached_indices and\
                [xind, yind + 1] not in self.solar_indices:
                    attached_indices.append([xind, yind + 1])
            if yind - 1 >= 0 and [xind, yind - 1] not in attached_indices and\
                [xind, yind - 1] not in self.solar_indices:
                    attached_indices.append([xind, yind - 1])
            if xind + 1 <= self.ncols - 1 and [xind + 1, yind] not in attached_indices and\
                [xind + 1, yind] not in self.solar_indices:
                    attached_indices.append([xind + 1, yind])
            if xind - 1 >= 0 and [xind - 1, yind] not in attached_indices and\
                [xind - 1, yind] not in self.solar_indices:
                    attached_indices.append([xind - 1, yind])
        max_distance = 0.0
        for i in range(len(attached_indices)):
            xind, yind = attached_indices[i]
            if self.distance[xind][yind] > max_distance:
                max_distance = self.distance[xind][yind]
                cell_index = [xind,yind]
                cell_location = self.grid_locs[xind][yind]
        self.distance[cell_index[0]][cell_index[1]] = -1E6
        
        if max_distance < self.min_spacing:
            return None, None
        else:
            return cell_location, cell_index
        

    def place_solar(self):
        self.reset()
        if self.nturbs > 0:
            self.find_distances()
            for i in range(self.nsolar_cells):
                if np.max(self.distance) >= self.min_spacing:
                    if type(self.solar_geometry) == geom.Polygon:
                        self.solar_geometry = geom.MultiPolygon([self.solar_geometry])
                    if len(self.solar_geometry) < self.max_arrays:
                        cell_location, cell_index = self.find_unattached_cell()
                    else:
                        cell_location, cell_index = self.find_attached_cell()

                    if cell_location == None:
                        break
                    
                    cell_geometry = geom.Polygon(([cell_location[0]-self.dx/2.0,cell_location[1]-self.dy/2.0],
                                                [cell_location[0]+self.dx/2.0,cell_location[1]-self.dy/2.0],
                                                [cell_location[0]+self.dx/2.0,cell_location[1]+self.dy/2.0],
                                                [cell_location[0]-self.dx/2.0,cell_location[1]+self.dy/2.0]))
                    if self.solar_geometry == None:
                        self.solar_geometry = cell_geometry
                    else:
                        self.solar_geometry = self.solar_geometry.union(cell_geometry.buffer(1E-12))
                    
                    self.solar_locs.append(cell_location)
                    self.solar_indices.append(cell_index)
                else:
                    break
        else:
            ncells = np.shape(self.grid_locs)[0] * np.shape(self.grid_locs)[1]
            if self.nsolar_cells > ncells:
                self.nsolar_cells = ncells
            
            total_area = self.dx*self.dy*self.nsolar_cells
            # print("total area: ", total_area)
            # print("self.nsolar_cells: ", self.nsolar_cells)
            side_length = np.sqrt(total_area)
            self.solar_geometry = geom.Polygon(([-side_length/2.0,-side_length/2.0],
                                                [side_length/2.0,-side_length/2.0],
                                                [side_length/2.0,side_length/2.0],
                                                [-side_length/2.0,side_length/2.0]))




if __name__=="__main__":

    high = 25
    x = np.linspace(0,1,high)
    y = np.linspace(0,1,high)
    # xlocs, ylocs = np.meshgrid(x,y)
    xlocs = [[i for i in x] for j in y]
    ylocs = [[j for i in x] for j in x]
    grid_locs = np.zeros((np.shape(xlocs)[0],np.shape(xlocs)[1],2))
    grid_locs[:,:,0] = xlocs[:]
    grid_locs[:,:,1] = ylocs[:]
    grid_locs = np.ndarray.tolist(grid_locs)
    print(np.shape(grid_locs))

    # turbine_x = [0.5]
    # turbine_y = [0.5]
    turbine_x = []
    turbine_y = []
    turbine_locs = np.zeros((len(turbine_x),2))
    turbine_locs[:,0] = turbine_x
    turbine_locs[:,1] = turbine_y
   

    nsolar_cells = 400
    max_arrays = 100
    min_spacing = 0.05

    start_solar = time.time()
    solar = PlaceSolar(grid_locs, min_spacing, max_arrays=max_arrays)
    solar.turbine_locs = turbine_locs
    solar.nsolar_cells = nsolar_cells
    solar.set_turbine_locs(turbine_locs)
    solar.place_solar()
    time_solar = time.time()-start_solar

    import matplotlib.pyplot as plt
    from plotting_functions import plot_poly, plot_turbines

    plot_turbines(turbine_x, turbine_y, 0.01, ax=plt.gca())
    plot_poly(solar.solar_geometry,ax=plt.gca())
    

    plt.axis("equal")
    plt.show()