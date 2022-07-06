# implement a simple shading model

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon, MultiPolygon
import time

class SimpleShadow():


    def __init__(self, hub_height, rotor_diameter, tower_width, latitide,  min_time=7, max_time=17):
        
        # turbine parameters
        self.hub_height = hub_height
        self.rotor_diameter = rotor_diameter
        self.tower_width = tower_width

        # position
        self.latitude = latitide

        self.min_time = min_time
        self.max_time = max_time

        self.turbine_x = np.array([])
        self.turbine_y = np.array([])
        self.solar_polygons = MultiPolygon()


    def rotate(self, origin, point, angle):
        """
        Rotate a point counterclockwise by a given angle around a given origin.

        The angle should be given in radians.
        """
        ox, oy = origin
        px, py = point

        qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (py - oy)
        qy = oy + np.sin(angle) * (px - ox) + np.cos(angle) * (py - oy)
        return qx, qy


    def find_angle(self, T_in):

        # find the omega
        from scipy import interpolate

        if T_in < 6:
            omega_out = 0
            print('Sun is not high enough for a shadow...')
        elif T_in > 18:
            omega_out = 0
            print('Sun is not high enough for a shadow...')
        else:
            omega_out = 15*T_in - 180.0

        return np.radians(omega_out)


    def calculate_shadow(self, t, day):

        # user inputs
        T = t  # time (in military time)
        d = day # number of days since the new year

        # turbine parameters
        HH = self.hub_height # hub height
        D = self.rotor_diameter # rotor diameter
        wd = self.tower_width  # tower width is 5 m?

        # position
        phi = np.radians(self.latitude)

        # calculate the shadow
        delta = np.radians(-23.45 * np.cos( np.radians(360/365 * (d + 10)) ))
        omega = self.find_angle(T)

 
        Fx = np.cos(delta) * np.sin(omega) / (np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.cos(omega))
        if Fx == 0.0:
            Fx = 1E-12
        numY = np.sin(phi) * np.cos(delta) * np.cos(omega) - np.cos(phi) * np.sin(delta)
        denY = np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.cos(omega) 
        Fy = numY / denY

        
        length = (HH + D/2) * Fx - (HH - D/2) * Fx
        angle = np.arctan2(Fy,Fx)

        a = length/2
        b = D/2
        N = 25
        x = np.linspace(-a,a,N)
        y = b * np.sqrt( 1 - (x/a)**2 )

        nturbs = len(self.turbine_x)
        tower_polys = np.zeros(nturbs,dtype=Polygon)
        rotor_polys = np.zeros(nturbs,dtype=Polygon)

        for k in range(nturbs):
            # turbine location
            x_loc = self.turbine_x[k]
            y_loc = self.turbine_y[k]

            poly_rotor = np.zeros((2*N,2))
            for i in range(len(x)):
                rx, ry = self.rotate([0,0], [x[i],y[i]], angle)
                # poly_rotor.append((rx[i]+(HH*Fx)+x_loc,ry[i]+(HH*Fy)+y_loc))
                poly_rotor[i][0] = rx+(HH*Fx)+x_loc
                poly_rotor[i][1] = ry+(HH*Fy)+y_loc
                # if t == 12.0:
                #     print("ry: ", ry)
                #     print(HH*Fy)

            for i in range(len(x)):
                rx2, ry2 = self.rotate([0, 0], [x[N-1-i], -y[N-1-i]], angle)
                poly_rotor[i+N][0] = rx2+(HH*Fx)+x_loc
                poly_rotor[i+N][1] = ry2+(HH*Fy)+y_loc
           
            poly_tower = [(x_loc + wd/2, y_loc), (x_loc - wd/2, y_loc),
                        (x_loc - wd/2 + HH * Fx, y_loc + HH * Fy), (x_loc + wd/2 + HH * Fx, y_loc + HH * Fy)]

            tower_polys[k] = Polygon(poly_tower)
            rotor_polys[k] = Polygon(poly_rotor)

        return tower_polys, rotor_polys


    def calculate_overlap(self, t, day):

        solar_area = self.solar_polygons.area
        
        rotor_ratio = 1./3.
        # rotor_ratio = 1.0

        tower_polys, rotor_polys = self.calculate_shadow(t,day)
            
        tower_shadows = tower_polys[0]
        for i in range(len(tower_polys)-1):
            tower_shadows = tower_shadows.union(tower_polys[i+1])
        
        tower_overlap = tower_shadows.buffer(0).intersection(self.solar_polygons.buffer(0))
        tower_loss = tower_overlap.area/solar_area

        rotor_shadows = rotor_polys[0]
        for i in range(len(rotor_polys)-1):
            rotor_shadows = rotor_shadows.union(rotor_polys[i+1])
        rotor_shadows = rotor_shadows.difference(tower_shadows)

        rotor_overlap = rotor_shadows.buffer(0).intersection(self.solar_polygons.buffer(0))
        rotor_loss = rotor_overlap.area*rotor_ratio/solar_area

        return tower_loss + rotor_loss


    def calculate_losses(self, t, day):
        if t < self.min_time or t > self.max_time:
            losses = 0.0
        else:
            if len(self.turbine_x) > 0:
                losses = self.calculate_overlap(t, day)
            else:
                losses = 0.0

        return losses


if __name__=="__main__":

        S = SimpleShadow()
        N = 1
        S.turbine_y = np.linspace(0.0,1000.0,N)
        S.turbine_x = np.zeros(N) + 10050.0
        S.latitude = 40.966
        poly = Polygon(((10100,0),(10100,1000),(10200,1000),(10200,0)))
        S.solar_polygons = MultiPolygon([poly])
        # S.latitude = 10.0
        

        start = time.time()

        n = 100
        t = np.linspace(0,24,n)
        # d = np.linspace(0,365,n)
        L = np.zeros(n)
        s = 80
        days = np.linspace(s,s+20,10)
        start = time.time()
        for k in range(len(days)):
            for i in range(n):
                L[i] = S.calculate_losses(t[i],days[k])
            plt.plot(t,L)
        plt.show()
        
        print("time to run: ", time.time()-start)        