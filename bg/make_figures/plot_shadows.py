import numpy as np
from shapely.geometry import Polygon
import matplotlib.pyplot as plt


def rotate(origin, point, angle):
        """
        Rotate a point counterclockwise by a given angle around a given origin.

        The angle should be given in radians.
        """
        ox, oy = origin
        px, py = point

        qx = ox + np.cos(angle) * (px - ox) - np.sin(angle) * (py - oy)
        qy = oy + np.sin(angle) * (px - ox) + np.cos(angle) * (py - oy)
        return qx, qy


def find_angle(T_in):

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


def calculate_shadow(t):

        T = t  # time (in military time)
        # d = 180.0 # number of days since the new year
        d = 180.0


        # turbine parameters
        HH = 200.0
        D = 135.0
        wd = 6.0

        # position
        phi = np.radians(40.966)
        # phi = np.radians(60.0)

        # calculate the shadow
        delta = np.radians(-23.45 * np.cos( np.radians(360/365 * (d + 10)) ))
        omega = find_angle(T)

 
        Fx = np.cos(delta) * np.sin(omega) / (np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.cos(omega))
        if Fx == 0.0:
            Fx = 1E-12
        numY = np.sin(phi) * np.cos(delta) * np.cos(omega) - np.cos(phi) * np.cos(delta)
        denY = np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.cos(omega) 
        Fy = numY / denY

        length = (HH + D/2) * Fx - (HH - D/2) * Fx
        angle = np.arctan2(Fy,Fx)

        a = length/2
        b = D/2
        N = 100
        x = np.linspace(-a,a,N)
        y = b * np.sqrt( 1 - (x/a)**2 )

        turbine_x = [0.0]
        turbine_y = [0.0]
        nturbs = len(turbine_x)
        tower_polys = np.zeros(nturbs,dtype=Polygon)
        rotor_polys = np.zeros(nturbs,dtype=Polygon)

        for k in range(nturbs):
            # turbine location
            x_loc = turbine_x[k]
            y_loc = turbine_y[k]

            poly_rotor = np.zeros((2*N,2))
            for i in range(len(x)):
                rx, ry = rotate([0,0], [x[i],y[i]], angle)
                # poly_rotor.append((rx[i]+(HH*Fx)+x_loc,ry[i]+(HH*Fy)+y_loc))
                poly_rotor[i][0] = rx+(HH*Fx)+x_loc
                poly_rotor[i][1] = ry+(HH*Fy)+y_loc

            for i in range(len(x)):
                rx2, ry2 = rotate([0, 0], [x[N-1-i], -y[N-1-i]], angle)
                # poly_rotor.append((rx2[i]+(HH*Fx)+x_loc,ry2[i]+(HH*Fy)+y_loc))
                poly_rotor[i+N][0] = rx2+(HH*Fx)+x_loc
                poly_rotor[i+N][1] = ry2+(HH*Fy)+y_loc
           

            poly_tower = [(x_loc + wd/2, y_loc), (x_loc - wd/2, y_loc),
                        (x_loc - wd/2 + (HH) * Fx, y_loc + HH * Fy), (x_loc + wd/2 + HH * Fx, y_loc + HH * Fy)]

            tower_polys[k] = Polygon(poly_tower)
            rotor_polys[k] = Polygon(poly_rotor)

        return tower_polys, rotor_polys

fig = plt.figure(figsize=[5,3])

t = 8.0
tower_polys, rotor_polys = calculate_shadow(t)
tx, ty = tower_polys[0].exterior.coords.xy
rx, ry = rotor_polys[0].exterior.coords.xy
plt.fill(tx,ty,color="C0")
plt.fill(rx,ry,color="C0",edgecolor=None,alpha=1./3.)
plt.text(rotor_polys[0].centroid.x,rotor_polys[0].centroid.y-30,"8:00",fontsize=8,horizontalalignment="center")

t = 10.0
tower_polys, rotor_polys = calculate_shadow(t)
tx, ty = tower_polys[0].exterior.coords.xy
rx, ry = rotor_polys[0].exterior.coords.xy
plt.fill(tx,ty,color="C0")
plt.fill(rx,ry,color="C0",edgecolor=None,alpha=1./3.)
plt.text(rotor_polys[0].centroid.x,rotor_polys[0].centroid.y+30,"10:00",fontsize=8,horizontalalignment="center")

t = 12.1
tower_polys, rotor_polys = calculate_shadow(t)
tx, ty = tower_polys[0].exterior.coords.xy
rx, ry = rotor_polys[0].exterior.coords.xy
plt.fill(tx,ty,color="C0")
plt.fill(rx,ry,color="C0",edgecolor=None,alpha=1./3.)
plt.text(rotor_polys[0].centroid.x,rotor_polys[0].centroid.y+30,"12:06",fontsize=8,horizontalalignment="center")

t = 14.0
tower_polys, rotor_polys = calculate_shadow(t)
tx, ty = tower_polys[0].exterior.coords.xy
rx, ry = rotor_polys[0].exterior.coords.xy
plt.fill(tx,ty,color="C0")
plt.fill(rx,ry,color="C0",edgecolor=None,alpha=1./3.)
plt.text(rotor_polys[0].centroid.x,rotor_polys[0].centroid.y+30,"14:00",fontsize=8,horizontalalignment="center")

t = 16.0
tower_polys, rotor_polys = calculate_shadow(t)
tx, ty = tower_polys[0].exterior.coords.xy
rx, ry = rotor_polys[0].exterior.coords.xy
plt.fill(tx,ty,color="C0")
plt.fill(rx,ry,color="C0",edgecolor=None,alpha=1./3.)
plt.text(rotor_polys[0].centroid.x,rotor_polys[0].centroid.y-30,"16:00",fontsize=8,horizontalalignment="center")

plt.axis("equal")
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['top'].set_visible(False)

plt.xlabel("x (m)", fontsize=8)
plt.ylabel("y (m)", fontsize=8)

plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

plt.title("shading model", fontsize=8)

plt.tight_layout()

# plt.savefig("figures/shading.pdf",transparent=True)
plt.show()