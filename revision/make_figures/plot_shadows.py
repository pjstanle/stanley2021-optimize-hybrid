import numpy as np
import matplotlib.pyplot as plt
from simple_shadow import SimpleShadow
from shapely.geometry import Polygon
from plotting_functions import plot_poly


turbine_x = np.array([0])
turbine_y = np.array([0])
solar_polygons = Polygon()
latitude = 40.966
hub_height = 90.0
rotor_diameter = 126.0
tower_width = 5.0

shadows = SimpleShadow(hub_height, rotor_diameter, tower_width, latitude,  min_time=7, max_time=17)
shadows.turbine_x = turbine_x
shadows.turbine_y = turbine_y

t = np.linspace(8,16,5)
day = 0

plt.figure(figsize=(3.75,5))
ax1 = plt.subplot(311)
ax2 = plt.subplot(312,sharex=ax1)
ax3 = plt.subplot(313,sharex=ax1)

ax1.set_title("January",fontsize=8)
ax2.set_title("February",fontsize=8)
ax3.set_title("March",fontsize=8)

for i in range(len(t)):
    tower_polys, rotor_polys = shadows.calculate_shadow(t[i], day)
    plot_poly(tower_polys[0],ax=ax1,color="C0",alpha=1.0)
    plot_poly(rotor_polys[0],ax=ax1,color="C0",alpha=1./3.)

ax1.axis("equal")

ax1.text(-500,550,"8:00",fontsize=6,horizontalalignment="center",verticalalignment="center")
ax1.text(-120,320,"10:00",fontsize=6,horizontalalignment="center",verticalalignment="center")
ax1.text(0,200,"12:00",fontsize=6,horizontalalignment="center",verticalalignment="center")
ax1.text(120,320,"14:00",fontsize=6,horizontalalignment="center",verticalalignment="center")
ax1.text(500,550,"16:00",fontsize=6,horizontalalignment="center",verticalalignment="center")

ax1.tick_params(axis='both', which='minor', labelsize=8)
ax1.tick_params(axis='both', which='major', labelsize=8)

ax1.set_ylim(-100,1010)

D = 30

for i in range(len(t)):
    tower_polys, rotor_polys = shadows.calculate_shadow(t[i], day+D)
    plot_poly(tower_polys[0],ax=ax2,alpha=1.0,color="C0")
    plot_poly(rotor_polys[0],ax=ax2,alpha=1./3.,color="C0")

ax2.axis("equal")

ax2.text(-500,550,"8:00",fontsize=6,horizontalalignment="center",verticalalignment="center")
ax2.text(-120,320,"10:00",fontsize=6,horizontalalignment="center",verticalalignment="center")
ax2.text(0,200,"12:00",fontsize=6,horizontalalignment="center",verticalalignment="center")
ax2.text(120,320,"14:00",fontsize=6,horizontalalignment="center",verticalalignment="center")
ax2.text(500,550,"16:00",fontsize=6,horizontalalignment="center",verticalalignment="center")

ax2.tick_params(axis='both', which='minor', labelsize=8)
ax2.tick_params(axis='both', which='major', labelsize=8)

ax2.set_ylim(-100,1010)

for i in range(len(t)):
    tower_polys, rotor_polys = shadows.calculate_shadow(t[i], day+2*D)
    plot_poly(tower_polys[0],ax=ax3,alpha=1.0,color="C0")
    plot_poly(rotor_polys[0],ax=ax3,alpha=1./3.,color="C0")

ax3.axis("equal")



ax3.text(-500,300,"8:00",fontsize=6,horizontalalignment="center",verticalalignment="center")
ax3.text(-120,200,"10:00",fontsize=6,horizontalalignment="center",verticalalignment="center")
ax3.text(0,400,"12:00",fontsize=6,horizontalalignment="center",verticalalignment="center")
ax3.text(120,200,"14:00",fontsize=6,horizontalalignment="center",verticalalignment="center")
ax3.text(500,300,"16:00",fontsize=6,horizontalalignment="center",verticalalignment="center")
arrow = plt.arrow(0,350,0,-250,)
ax3.add_patch(arrow)

ax3.tick_params(axis='both', which='minor', labelsize=8)
ax3.tick_params(axis='both', which='major', labelsize=8)

ax3.set_ylim(-100,1010)

ax1.set_yticks((0,500,1000))
ax2.set_yticks((0,500,1000))
ax3.set_yticks((0,500,1000))

ax3.set_xlabel("x coordinate",fontsize=8)
ax1.set_ylabel("y coordinate",fontsize=8)
ax2.set_ylabel("y coordinate",fontsize=8)
ax3.set_ylabel("y coordinate",fontsize=8)

# plt.tight_layout()
plt.subplots_adjust(top=0.95,bottom=0.1,left=0.2,right=0.95,hspace=0.6)
plt.savefig("shadows_mat.pdf",transparent=True)
plt.show()
