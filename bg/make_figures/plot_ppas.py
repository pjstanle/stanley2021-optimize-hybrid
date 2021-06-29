import numpy as np
import matplotlib.pyplot as plt
from boundary_grid import makeBoundary, makeGrid
from shapely.geometry import Polygon, Point, LineString, MultiPolygon


def get_solar_polygon(solar_x,solar_y,solar_width,solar_height):

    solar_pts = [(solar_x-solar_width/2.0,solar_y-solar_height/2.0),
                (solar_x+solar_width/2.0,solar_y-solar_height/2.0),
                (solar_x+solar_width/2.0,solar_y+solar_height/2.0),
                (solar_x-solar_width/2.0,solar_y+solar_height/2.0),
                (solar_x-solar_width/2.0,solar_y-solar_height/2.0)]
    
    pts_x = np.array([solar_x-solar_width/2.0,solar_x+solar_width/2.0,solar_x+solar_width/2.0,solar_x-solar_width/2.0])
    pts_y = np.array([solar_y-solar_height/2.0,solar_y-solar_height/2.0,solar_y+solar_height/2.0,solar_y+solar_height/2.0])
    
    solar_poly = Polygon(solar_pts)
    return solar_poly, pts_x, pts_y


def combine_solar_polygons(poly0,poly1,poly2):
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

    return solar_union


def delete_solar_turbines(initial_turbine_x,initial_turbine_y,solar_union):
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
    
    return turbine_x, turbine_y


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


def get_xy(A):
    x = np.zeros(len(A))
    y = np.zeros(len(A))
    for i in range(len(A)):
        x[i] = A[i][0]
        y[i] = A[i][1]
    return x,y


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
    # try:
    #     for j in range(len(solar_polys)):
    #         x,y = solar_polys[j].boundary.coords.xy
    #         plt.fill(x,y,color="C1")
    # except:
    #     x,y = solar_polys.boundary.coords.xy
    #     plt.fill(x,y,color="C1")
    if type(solar_polys) == Polygon:
        x,y = solar_polys.boundary.coords.xy
        plt.fill(x,y,color="C1")
    else:
        for part in solar_polys:
            exterior_coords = part.exterior.coords[:]
            x,y = get_xy(exterior_coords)
            plt.fill(x,y,color="C1")

            for interior in part.interiors:
                interior_coords = interior.coords[:]
                x,y = get_xy(interior_coords)
                plt.fill(x,y,"white")


    # plot turbines
    for i in range(len(turbine_x)):
        # turb = plt.Circle((turbine_x[i],turbine_y[i]),rotor_diameter/2.0,facecolor="C0")
        turb = plt.Circle((turbine_x[i],turbine_y[i]),rotor_diameter/2.0,facecolor="C0")
        plt.gca().add_patch(turb)


def get_dvs(inputs):

    boundary_start = inputs[0]
    boundary_turbs = inputs[1]
    grid_rows = inputs[2]
    grid_cols = inputs[3]
    grid_spacing_x = inputs[4]
    grid_spacing_y = inputs[5]
    grid_shear = inputs[6]
    grid_rotation = inputs[7]
    grid_center_x = inputs[8]
    grid_center_y = inputs[9]
    grid_boundary_setback = inputs[10]

    solar_x = inputs[11]
    solar_y = inputs[12]
    solar_width = inputs[13]
    solar_height = inputs[14]

    # battery variables
    battery_storage = inputs[15]

    # assign turbine variables
    boundary_turbine_x, boundary_turbine_y = makeBoundary(boundary_start,boundary_turbs,boundary_poly)
    grid_turbine_x, grid_turbine_y = makeGrid(grid_rows,grid_cols,grid_spacing_x,grid_spacing_y,grid_shear,
                                                grid_rotation,grid_center_x,grid_center_y,grid_boundary_setback,boundary_poly)
    initial_turbine_x = np.append(boundary_turbine_x,grid_turbine_x)
    initial_turbine_y = np.append(boundary_turbine_y,grid_turbine_y)

    poly,pts_x,pts_y = get_solar_polygon(solar_x,solar_y,solar_width,solar_height)
    poly = poly.intersection(boundary_poly)

   
    if type(poly)==Polygon:
        solar_union = MultiPolygon([poly])
    else:
        solar_union = poly

    turbine_x, turbine_y = delete_solar_turbines(initial_turbine_x,initial_turbine_y,solar_union)

    return turbine_x, turbine_y, solar_union


global boundary_poly
global rotor_diameter
global wind_solar_spacing
global wind_wind_spacing

rotor_diameter = 200.0
D = rotor_diameter
wind_wind_spacing = D*4.0
wind_solar_spacing = D*0.5
boundary_poly = Polygon(([0,0],[1000,-1000],[1200,1000],[500,1100],[-200,1000])).buffer(2000)

time = np.linspace(6,30,100)

plt.figure(figsize=(5,5.5))

plt.subplot(331)
# PPA
def PPA_curve1(hour):
    return 70.0
ppa = np.zeros(len(time))
for i in range(len(time)):
    ppa[i] = PPA_curve1(time[i])
plt.plot(time,ppa)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.ylim(59,71)
plt.ylabel("PPA ($/MWh)",fontsize=8,labelpad=0)
plt.gca().set_yticks((60,65,70))
plt.gca().set_yticklabels(("60","65","70"))
plt.gca().set_xticks((6,13,24,30))
plt.gca().set_xticklabels(("6:00","13:00","24:00","6:00"))
plt.xlabel("time of day",labelpad=0,fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)

plt.subplot(332)
def PPA_curve2(hour):
    time = hour%24
    p1 = np.cos((13-time)/(2*np.pi))*70.0
    p2 = 60.0
    return max([p1,p2])
ppa = np.zeros(len(time))
for i in range(len(time)):
    ppa[i] = PPA_curve2(time[i])
plt.plot(time,ppa)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.ylim(59,71)
plt.gca().set_yticks((60,65,70))
plt.gca().set_yticklabels(("","",""))
plt.gca().set_xticks((6,13,24,30))
plt.gca().set_xticklabels(("6:00","13:00","24:00","6:00"))
plt.xlabel("time of day",labelpad=0,fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
print(time)
print(ppa)


plt.subplot(333)
def PPA_curve3(hour):
    time = hour%24
    p1 = np.cos((time)/(2*np.pi))*70.0
    p2 = np.cos((24-time)/(2*np.pi))*70.0
    p3 = 60.0
    return max([p1,p2,p3])
ppa = np.zeros(len(time))
for i in range(len(time)):
    ppa[i] = PPA_curve3(time[i])
plt.plot(time,ppa)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.ylim(59,71)
plt.gca().set_yticks((60,65,70))
plt.gca().set_yticklabels(("","",""))
plt.gca().set_xticks((6,13,24,30))
plt.gca().set_xticklabels(("6:00","13:00","24:00","6:00"))
plt.xlabel("time of day",labelpad=0,fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)



plt.subplot(334)
# 1
# solar capacity:  309.96486946935073
# wind capacity:  231.0
# battery capacity:  361.9047619047619
# coe:  -17.60228961599387

plt.bar([1],[231],width=0.5,color="C0",bottom=309.96486946935073,label="wind (MW)")
plt.bar([1],[309.96486946935073],width=0.5,color="C1",label="solar (MW)")
plt.plot([4],[361.9047619047619],"o",color="C3",label="battery (MWh)")
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlim(0,5)
plt.ylim(0,599)
plt.legend(loc=4,fontsize=6)
plt.gca().set_xticks((1,4))
plt.gca().set_xticklabels(("generation","storage"),fontsize=8)
plt.gca().set_yticks((0,250,500))
plt.gca().set_yticklabels(("0","250","500"),fontsize=8)
plt.ylabel("capacity",fontsize=8,labelpad=0)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)


plt.subplot(335)
# 2
# solar capacity:  316.45178728640764
# wind capacity:  217.0
# battery capacity:  361.9047619047619
# coe:  -6.398059000381083

plt.bar([1],[316.45178728640764],width=0.5,color="C1")
plt.bar([1],[217],width=0.5,color="C0",bottom=316.45178728640764)
plt.plot([4],[361.9047619047619],"o",color="C3")
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlim(0,5)
plt.ylim(0,599)
plt.gca().set_xticks((1,4))
plt.gca().set_xticklabels(("generation","storage"),fontsize=8)
plt.gca().set_yticks((0,250,500))
plt.gca().set_yticklabels(("","",""),fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)


plt.subplot(336)
# 3
# solar capacity:  295.2046375898579
# wind capacity:  252.0
# battery capacity:  361.9047619047619
# coe:  -3.9609444021381437

plt.bar([1],[295.2046375898579],width=0.5,color="C1")
plt.bar([1],[252],width=0.5,color="C0",bottom=295.2046375898579)
plt.plot([4],[361.9047619047619],"o",color="C3")
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlim(0,5)
plt.ylim(0,599)
plt.gca().set_xticks((1,4))
plt.gca().set_xticklabels(("generation","storage"),fontsize=8)
plt.gca().set_yticks((0,250,500))
plt.gca().set_yticklabels(("","",""),fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)



plt.subplot(337)
inputs = np.array([ 1.1746031746031747e+03,  2.2000000000000000e+01,
        1.2000000000000000e+01,  1.5000000000000000e+01,
        8.1587301587301590e+02,  1.0317460317460318e+03,
       -8.2279807594018384e-01,  4.6874557053561992e+00,
        1.5093319986198167e+03, -1.3853167398297160e+03,
        7.6190476190476193e+02,  1.3791711312916550e+02,
       -3.0258306510178500e+01,  3.0856834923539654e+03,
        2.0325876499793058e+03,  3.6190476190476193e+02])
d1,d2,d3 = get_dvs(inputs)
print(len(d1)*7)
plot_hybrid(d1,d2,d3)
plt.axis("equal")
plt.axis("off")
plt.text(np.mean(d1),np.min(d2)-200,"$17.60 MM",fontsize=8,horizontalalignment="center",verticalalignment="top")

plt.subplot(338)
# 2
inputs = np.array([ 6.6666666666666674e+02,  2.2000000000000000e+01,
    1.2000000000000000e+01,  1.3000000000000000e+01,
    8.1587301587301590e+02,  7.0793650793650795e+02,
    5.2359877559829915e-01,  4.9866550056980845e+00,
    2.1093260110219767e+03, -8.0457741126419978e+02,
    8.7301587301587301e+02,  1.0807648469039882e+03,
    2.8734383363174020e+03,  2.0571223282359770e+03,
    5.9525781177965382e+03,  3.6190476190476193e+02])
d1,d2,d3 = get_dvs(inputs)
print(len(d1)*7)
plot_hybrid(d1,d2,d3)
plt.axis("equal")
plt.axis("off")
plt.text(np.mean(d1),np.min(d2)-200,"$6.40 MM",fontsize=8,horizontalalignment="center",verticalalignment="top")

plt.subplot(339)
# 3
inputs = np.array([ 6.0317460317460313e+02,  2.2000000000000000e+01,
    1.1000000000000000e+01,  1.2000000000000000e+01,
    8.1587301587301590e+02,  8.1587301587301590e+02,
   -2.7426602531339461e-01,  1.1967972013675403e+00,
    1.5950454289629824e+03, -9.0136729935845256e+02,
    7.6190476190476193e+02,  6.5219769518815974e+02,
   -3.0258306510178500e+01,  4.1142446564719539e+03,
    1.4518483214137900e+03,  3.6190476190476193e+02])
d1,d2,d3 = get_dvs(inputs)
print(len(d1)*7)
plot_hybrid(d1,d2,d3)
plt.axis("equal")
plt.axis("off")
plt.text(np.mean(d1),np.min(d2)-200,"$3.96 MM",fontsize=8,horizontalalignment="center",verticalalignment="top")


plt.subplots_adjust(top=0.99,bottom=0.02,left=0.1,right=0.98,hspace=0.3)

# plt.savefig("figures/ppa_table.pdf",transparent=True)
plt.show()