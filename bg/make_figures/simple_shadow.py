# implement a simple shading model

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon, MultiPolygon
import time

class SimpleShadow():


    def __init__(self):

        self.turbine_x = np.array([])
        self.turbine_y = np.array([])
        self.solar_x = np.array([]) # coordinates of the solar array vertices. Right now assuming 4 vertices
        self.solar_y = np.array([])
        
        # turbine parameters
        self.hub_height = 90.0
        self.rotor_diameter = 126.4
        self.tower_width = 5.0

        # position
        self.latitude = 39.7555
        self.longitude = -105.2211


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


    def calculate_shadow(self, t, day, rotor_direction):

        global delta_array
        global omega_array
        global Fx_array
        global numY_array
        global denY_array
        global Fy_array

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
        numY = np.sin(phi) * np.cos(delta) * np.cos(omega) - np.cos(phi) * np.cos(delta)
        denY = np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.cos(omega) 
        Fy = numY / denY

        length = (HH + D/2) * Fx - (HH - D/2) * Fx
        angle = np.arctan2(Fy,Fx)

        # adjust the width of the shadow in case the turbine is not facing directly into the sun
        # this works fine (I think) except when the sun is close to directly overhead, in which case it shrinks the shadow more than it should
        # width_multiplier = abs(np.sin(angle-np.radians(rotor_direction))) 

        a = length/2
        # b = D/2 * width_multiplier
        b = D/2
        N = 25
        x = np.linspace(-a,a,N)
        y = b * np.sqrt( 1 - (x/a)**2 )

        # rx = np.zeros(len(x))
        # ry = np.zeros(len(y))
        # rx2 = np.zeros(len(x))
        # ry2 = np.zeros(len(y))

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

            for i in range(len(x)):
                rx2, ry2 = self.rotate([0, 0], [x[N-1-i], -y[N-1-i]], angle)
                # poly_rotor.append((rx2[i]+(HH*Fx)+x_loc,ry2[i]+(HH*Fy)+y_loc))
                poly_rotor[i+N][0] = rx2+(HH*Fx)+x_loc
                poly_rotor[i+N][1] = ry2+(HH*Fy)+y_loc
           

            poly_tower = [(x_loc + wd/2, y_loc), (x_loc - wd/2, y_loc),
                        (x_loc - wd/2 + (HH) * Fx, y_loc + HH * Fy), (x_loc + wd/2 + HH * Fx, y_loc + HH * Fy)]

            tower_polys[k] = Polygon(poly_tower)
            rotor_polys[k] = Polygon(poly_rotor)

        return tower_polys, rotor_polys


    def calculate_overlap(self, t, day, rotor_direction=0, show=False):

        # determine xmin, xmax, ymin, ymax
        # assume 4 sided

        s = time.time()
        solar_poly = Polygon(((self.solar_x[0],self.solar_y[0]),(self.solar_x[1],self.solar_y[1]),(self.solar_x[2],self.solar_y[2]),(self.solar_x[3],self.solar_y[3])))
        solar_area = solar_poly.area
        
        rotor_ratio = 1./3.
        # rotor_ratio = 1.0
        # print("create solar poly: ", time.time()-s)

        s = time.time()
        tower_polys, rotor_polys = self.calculate_shadow(t,day,rotor_direction)
        # print("create shadow polys: ", time.time()-s)
            
        s = time.time()
        tower_shadows = tower_polys[0]
        for i in range(len(tower_polys)-1):
            tower_shadows = tower_shadows.union(tower_polys[i+1])
        
        tower_overlap = tower_shadows.intersection(solar_poly)
        tower_loss = tower_overlap.area/solar_area
        # print("loss towers: ", time.time()-s)

        s = time.time()
        rotor_shadows = rotor_polys[0]
        for i in range(len(rotor_polys)-1):
            rotor_shadows = rotor_shadows.union(rotor_polys[i+1])
        rotor_shadows = rotor_shadows.difference(tower_shadows)

        rotor_overlap = rotor_shadows.intersection(solar_poly)
        rotor_loss = rotor_overlap.area*rotor_ratio/solar_area
        # print("union rotor: ", time.time()-s)


        if show:
            if type(tower_shadows) == Polygon:
                tower_shadows = MultiPolygon([tower_shadows])
            if type(rotor_shadows) == Polygon:
                rotor_shadows = MultiPolygon([rotor_shadows])
            if type(tower_overlap) == Polygon:
                tower_overlap = MultiPolygon([tower_overlap])
            if type(rotor_overlap) == Polygon:
                rotor_overlap = MultiPolygon([rotor_overlap])

            plt.gca().cla()
            x,y = solar_poly.exterior.coords.xy
            plt.fill(x,y,"b",label="solar")
            if tower_shadows.area > 0.0:
                for i in range(len(tower_shadows)):
                    x,y = tower_shadows[i].exterior.coords.xy
                    plt.fill(x,y,"k")
            if tower_overlap.area > 0.0:
                for i in range(len(tower_overlap)):
                    x,y = tower_overlap[i].exterior.coords.xy
                    plt.fill(x,y,"r")
            if rotor_shadows.area > 0.0:
                for i in range(len(rotor_shadows)):
                    x,y = rotor_shadows[i].exterior.coords.xy
                    if i==0:
                        plt.fill(x,y,"k",label="shadow")
                    else:
                        plt.fill(x,y,"k")
            if rotor_overlap.area > 0.0:
                for i in range(len(rotor_overlap)):
                    x,y = rotor_overlap[i].exterior.coords.xy
                    if i==0:
                        plt.fill(x,y,"r",label="overlap")
                    else:
                        plt.fill(x,y,"r")
            
            plt.legend()
            plt.axis("square")
            plt.ylim(-2000,2000)
            plt.xlim(8000,11000)
            plt.pause(0.2)

        return tower_loss + rotor_loss


    def calculate_losses(self, t, day, rotor_direction=0, show_shadow=False):
        if len(self.turbine_x) > 0:
            losses = self.calculate_overlap(t, day, rotor_direction=rotor_direction, show=show_shadow)
        else:
            losses = 0.0
        return losses


if __name__=="__main__":

        global delta_array
        global omega_array
        global Fx_array
        global numY_array
        global denY_array
        global Fy_array

        delta_array = np.array([])
        omega_array = np.array([])
        Fx_array = np.array([])
        numY_array = np.array([])
        denY_array = np.array([])
        Fy_array = np.array([])

        S = SimpleShadow()
        N = 25
        S.turbine_y = np.linspace(0.0,1000.0,N)
        S.turbine_x = np.zeros(N) + 10050.0
        S.solar_x = np.array([9000.0,10000.0,10000.0,9000.0])
        S.solar_y = np.array([0.0,0.0,1000.0,1000.0])
        S.latitude = 40.966
        # S.latitude = 10.0
        

        start = time.time()

        n = 1
        t = np.linspace(7,17,n)
        d = np.linspace(0,365,n)
        L = np.zeros(n)
        start = time.time()
        for i in range(n):
            print(t[i])
            # L[i] = S.calculate_losses(13.0,d[i],rotor_direction=0.0,show_shadow=False)
            L[i] = S.calculate_losses(t[i],10.0,rotor_direction=0.0,show_shadow=False)
            # print("time to run: ", time.time()-start)
            # print("loss: ", L)
        
        print("time to run: ", time.time()-start)

        # plt.gca().cla()
        print(L)
        plt.plot(L)
        plt.show()