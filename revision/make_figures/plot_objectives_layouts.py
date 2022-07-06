import numpy as np
import matplotlib.pyplot as plt
from plotting_functions import plot_poly, plot_turbines
from make_solar import PlaceSolar
from shapely.geometry import Polygon


if __name__=="__main__":

    corner = 2600.0
    boundary_polygon = Polygon(((-corner,-corner),(-corner,corner),(corner,corner),(corner,-corner)))

    powercurve_filename = 'high_7r_200d_135h.txt'
    rotor_diameter = 200.0
    hub_height = 135.0
    turbine_rating = 7.0

    corner = 2500.0
    ngrid = 25 # dunno what this should be, affects how coarse the solar grid is
    x = np.linspace(-corner,corner,ngrid)
    y = np.linspace(-corner,corner,ngrid)
    xlocs = [[i for i in x] for j in y]
    ylocs = [[j for i in x] for j in x]
    grid_locs = np.zeros((np.shape(xlocs)[0],np.shape(xlocs)[1],2))
    grid_locs[:,:,0] = xlocs[:]
    grid_locs[:,:,1] = ylocs[:]
    grid_locs = np.ndarray.tolist(grid_locs)

    dx = grid_locs[0][1][0] - grid_locs[0][0][0]
    dy = grid_locs[1][0][1] - grid_locs[0][0][1]
    solar_cell_area = dx*dy
    solar_kw_per_km2 = 1000.0/5.0 * 247.105

    min_solar_spacing = 2*rotor_diameter
    place_solar = PlaceSolar(grid_locs,min_solar_spacing)



    plt.figure(figsize=(4,5.75))
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)

    plot_poly(boundary_polygon,ax=ax1,color="black",fill=False)
    plot_poly(boundary_polygon,ax=ax2,color="black",fill=False)
    plot_poly(boundary_polygon,ax=ax3,color="black",fill=False)
    plot_poly(boundary_polygon,ax=ax4,color="black",fill=False)


    # not_resil_const
    wind = 112000.
    solar = 0.0
    battery = 0.0

    coe = 34.54051255673027
    profit = 16271041.008985266
    capacity_factor = 0.1746053413478216

    turbine_x = np.array([ 1875.00676725,  2383.33333333,  -252.05205628,   256.2745098 ,
         764.60107589,  1272.92764197,  1781.25420806,  2289.58077414,
       -2379.11087981, -1870.78431373, -1362.45774764,  -854.13118156,
        -345.80461547,   162.52195061, -2472.863439  , -1964.53687292])

    turbine_y = np.array([ 2496.66666667,  1499.01960784,  2496.66666667,  1499.01960784,
         501.37254902,  -496.2745098 , -1493.92156863, -2491.56862745,
        2496.66666667,  1499.01960784,   501.37254902,  -496.2745098 ,
       -1493.92156863, -2491.56862745, -1493.92156863, -2491.56862745])

          
    nturbs = len(turbine_x)
    turbine_locs = np.zeros((nturbs,2))
    turbine_locs[:,0] = turbine_x
    turbine_locs[:,1] = turbine_y
    place_solar.set_turbine_locs(turbine_locs)
    place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
    place_solar.place_solar()
    solar_poly = place_solar.solar_geometry

    plot_poly(solar_poly,ax=ax1,color="C1")
    plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax1)

    DY = 3200
    ax1.text(0,0-DY,"COE: %s"%round(coe,2),fontsize=8,horizontalalignment="center")
    ax1.text(0,-500-DY,"profit: %s"%round(profit/1E6,2),fontsize=8,horizontalalignment="center")
    ax1.text(0,-1000-DY,"solar: %s"%round(solar/1E3,1),fontsize=8,horizontalalignment="center")
    ax1.text(0,-1500-DY,"wind: %s"%round(wind/1E3,1),fontsize=8,horizontalalignment="center")
    ax1.text(0,-2000-DY,"battery: %s"%round(battery/1E3,1),fontsize=8,horizontalalignment="center")



    # not_resil_const
    profit = 37077516.13549378
    wind = 301000.0
    solar = 42900.173611111626
    battery = 0.0

    coe = 36.92644114022708
    profit = 37077516.39740612
    capacity_factor = 0.42658378696574484

    turbine_x = np.array([ 2427.92279693,  2105.47794305,  1599.68750281,  2288.82352941,
            1783.03308917,  1277.24264893,   771.45220869,  2472.16911577,
            1966.37867553,  1460.58823529,   954.79779506,   449.00735482,
            -56.78308542,  2149.72426189,  1643.93382165,  1138.14338142,
            632.35294118,   126.56250094,  -379.2279393 ,  -885.01837954,
            1321.48896778,   815.69852754,   309.9080873 ,  -195.88235294,
            -701.67279318, -1207.46323342, -1713.25367366,   493.25367366,
            -12.53676658,  -518.32720682, -1024.11764706, -1529.9080873 ,
        -2035.69852754,  -334.98162046,  -840.7720607 , -1346.56250094,
        -1852.35294118, -2358.14338142, -1163.21691458, -1669.00735482,
        -2174.79779506, -1991.45220869, -2497.24264893])

    turbine_y = np.array([-2470.        , -1641.76470588, -2470.        ,    14.70588235,
            -813.52941176, -1641.76470588, -2470.        ,  1671.17647059,
            842.94117647,    14.70588235,  -813.52941176, -1641.76470588,
        -2470.        ,  2499.41176471,  1671.17647059,   842.94117647,
            14.70588235,  -813.52941176, -1641.76470588, -2470.        ,
            2499.41176471,  1671.17647059,   842.94117647,    14.70588235,
            -813.52941176, -1641.76470588, -2470.        ,  2499.41176471,
            1671.17647059,   842.94117647,    14.70588235,  -813.52941176,
        -1641.76470588,  2499.41176471,  1671.17647059,   842.94117647,
            14.70588235,  -813.52941176,  2499.41176471,  1671.17647059,
            842.94117647,  2499.41176471,  1671.17647059])

          
    nturbs = len(turbine_x)
    turbine_locs = np.zeros((nturbs,2))
    turbine_locs[:,0] = turbine_x
    turbine_locs[:,1] = turbine_y
    place_solar.set_turbine_locs(turbine_locs)
    place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
    place_solar.place_solar()
    solar_poly = place_solar.solar_geometry

    plot_poly(solar_poly,ax=ax2,color="C1")
    plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax2)

    DY = 3200
    ax2.text(0,0-DY,"COE: %s"%round(coe,2),fontsize=8,horizontalalignment="center")
    ax2.text(0,-500-DY,"profit: %s"%round(profit/1E6,2),fontsize=8,horizontalalignment="center")
    ax2.text(0,-1000-DY,"solar: %s"%round(solar/1E3,1),fontsize=8,horizontalalignment="center")
    ax2.text(0,-1500-DY,"wind: %s"%round(wind/1E3,1),fontsize=8,horizontalalignment="center")
    ax2.text(0,-2000-DY,"battery: %s"%round(battery/1E3,1),fontsize=8,horizontalalignment="center")



    # resil coe
    # wind = 224000.
    # solar = 349636.41493056
    # battery = 218963.83186706

    # coe = 63.88388786254236
    # profit = 5535772.651028618

    # turbine_x = np.array([-2.45741609e+03, -1.23393311e+03, -2.45231421e+03, -1.22883122e+03,
    #    -5.34823684e+00,  1.21813475e+03,  2.44161773e+03, -2.44721232e+03,
    #    -1.22372934e+03, -2.46351625e-01,  1.22323663e+03,  2.44671962e+03,
    #    -2.44211044e+03, -1.21862745e+03,  4.85553359e+00,  1.22833852e+03,
    #     2.45182150e+03, -2.43700855e+03, -1.21352557e+03,  9.95741881e+00,
    #     1.23344040e+03,  2.45692339e+03, -2.43190667e+03, -1.20842368e+03,
    #     1.50593040e+01,  1.23854229e+03,  2.46202527e+03, -2.42680478e+03,
    #    -1.20332180e+03,  2.01611892e+01,  1.24364417e+03,  2.46712716e+03])

    # turbine_y = np.array([ 2.48413778e+03,  2.49544305e+03,  1.65591820e+03,  1.66722347e+03,
    #     1.67852875e+03,  1.68983403e+03,  1.70113931e+03,  8.27698616e+02,
    #     8.39003894e+02,  8.50309172e+02,  8.61614450e+02,  8.72919728e+02,
    #    -5.20964356e-01,  1.07843137e+01,  2.20895918e+01,  3.33948699e+01,
    #     4.47001480e+01, -8.28740545e+02, -8.17435267e+02, -8.06129988e+02,
    #    -7.94824710e+02, -7.83519432e+02, -1.65696012e+03, -1.64565485e+03,
    #    -1.63434957e+03, -1.62304429e+03, -1.61173901e+03, -2.48517971e+03,
    #    -2.47387443e+03, -2.46256915e+03, -2.45126387e+03, -2.43995859e+03])

    profit = 5645419.810419798
    wind = 231000.0
    solar = 296011.1979166696
    battery = 222873.90029325514

    coe = 63.86832137209415
    profit = 5645420.447551072
    capacity_factor = 0.35034146878732314

    turbine_x = np.array([-2326.83074483, -2326.83074483, -2326.83074483, -2326.83074483,
        -2326.83074483, -2326.83074483, -1132.35294118, -1132.35294118,
        -1132.35294118, -1132.35294118, -1132.35294118, -1132.35294118,
            62.12486248,    62.12486248,    62.12486248,    62.12486248,
            62.12486248,    62.12486248,    62.12486248,  1256.60266613,
            1256.60266613,  1256.60266613,  1256.60266613,  1256.60266613,
            1256.60266613,  1256.60266613,  2451.08046978,  2451.08046978,
            2451.08046978,  2451.08046978,  2451.08046978,  2451.08046978,
            2451.08046978])

    turbine_y = np.array([-1761.06069328,  -960.51396801,  -159.96724274,   640.57948252,
            1441.12620779,  2241.67293306, -1716.89278342,  -916.34605815,
            -115.79933289,   684.74739238,  1485.29411765,  2285.84084291,
        -2473.27159883, -1672.72487356,  -872.17814829,   -71.63142303,
            728.91530224,  1529.46202751,  2330.00875277, -2429.10368897,
        -1628.5569637 ,  -828.01023843,   -27.46351317,   773.0832121 ,
            1573.62993736,  2374.17666263, -2384.93577911, -1584.38905384,
            -783.84232858,    16.70439669,   817.25112196,  1617.79784722,
            2418.34457249])

    nturbs = len(turbine_x)
    turbine_locs = np.zeros((nturbs,2))
    turbine_locs[:,0] = turbine_x
    turbine_locs[:,1] = turbine_y
    place_solar.set_turbine_locs(turbine_locs)
    place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
    place_solar.place_solar()
    solar_poly = place_solar.solar_geometry

    plot_poly(solar_poly,ax=ax3,color="C1")
    plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax3)

    DY = 3200
    ax3.text(0,0-DY,"COE: %s"%round(coe,2),fontsize=8,horizontalalignment="center")
    ax3.text(0,-500-DY,"profit: %s"%round(profit/1E6,2),fontsize=8,horizontalalignment="center")
    ax3.text(0,-1000-DY,"solar: %s"%round(solar/1E3,1),fontsize=8,horizontalalignment="center")
    ax3.text(0,-1500-DY,"wind: %s"%round(wind/1E3,1),fontsize=8,horizontalalignment="center")
    ax3.text(0,-2000-DY,"battery: %s"%round(battery/1E3,1),fontsize=8,horizontalalignment="center")



    # resil_const
    profit = 5645419.810419798
    wind = 231000.0
    solar = 296011.1979166696
    battery = 222873.90029325514
    capacity_factor = 0.35034146878732314

    coe = 63.86832137209415
    profit = 5645420.447551072

    turbine_x = np.array([-2326.83074483, -2326.83074483, -2326.83074483, -2326.83074483,
        -2326.83074483, -2326.83074483, -1132.35294118, -1132.35294118,
        -1132.35294118, -1132.35294118, -1132.35294118, -1132.35294118,
            62.12486248,    62.12486248,    62.12486248,    62.12486248,
            62.12486248,    62.12486248,    62.12486248,  1256.60266613,
            1256.60266613,  1256.60266613,  1256.60266613,  1256.60266613,
            1256.60266613,  1256.60266613,  2451.08046978,  2451.08046978,
            2451.08046978,  2451.08046978,  2451.08046978,  2451.08046978,
            2451.08046978])

    turbine_y = np.array([-1761.06069328,  -960.51396801,  -159.96724274,   640.57948252,
            1441.12620779,  2241.67293306, -1716.89278342,  -916.34605815,
            -115.79933289,   684.74739238,  1485.29411765,  2285.84084291,
        -2473.27159883, -1672.72487356,  -872.17814829,   -71.63142303,
            728.91530224,  1529.46202751,  2330.00875277, -2429.10368897,
        -1628.5569637 ,  -828.01023843,   -27.46351317,   773.0832121 ,
            1573.62993736,  2374.17666263, -2384.93577911, -1584.38905384,
            -783.84232858,    16.70439669,   817.25112196,  1617.79784722,
            2418.34457249])

          
    nturbs = len(turbine_x)
    turbine_locs = np.zeros((nturbs,2))
    turbine_locs[:,0] = turbine_x
    turbine_locs[:,1] = turbine_y
    place_solar.set_turbine_locs(turbine_locs)
    place_solar.nsolar_cells = int(np.floor(solar/(solar_kw_per_km2*solar_cell_area/1E6)))
    place_solar.place_solar()
    solar_poly = place_solar.solar_geometry

    plot_poly(solar_poly,ax=ax4,color="C1")
    plot_turbines(turbine_x, turbine_y, rotor_diameter/2.0, ax=ax4)

    DY = 3200
    ax4.text(0,0-DY,"COE: %s"%round(coe,2),fontsize=8,horizontalalignment="center")
    ax4.text(0,-500-DY,"profit: %s"%round(profit/1E6,2),fontsize=8,horizontalalignment="center")
    ax4.text(0,-1000-DY,"solar: %s"%round(solar/1E3,1),fontsize=8,horizontalalignment="center")
    ax4.text(0,-1500-DY,"wind: %s"%round(wind/1E3,1),fontsize=8,horizontalalignment="center")
    ax4.text(0,-2000-DY,"battery: %s"%round(battery/1E3,1),fontsize=8,horizontalalignment="center")




    ax1.set_title("minimize COE",fontsize=8)
    ax2.set_title("maximize profit",fontsize=8)
    ax1.text(-3200,0,"no minimum\npower constraint",fontsize=8,rotation=90,
        horizontalalignment="center",verticalalignment="center")
    ax3.text(-3200,0,"10-MW minimum\npower constraint",fontsize=8,rotation=90,
        horizontalalignment="center",verticalalignment="center")

    ax1.axis("square")
    ax1.axis("off")
    ax2.axis("square")
    ax2.axis("off")
    ax3.axis("square")
    ax3.axis("off")
    ax4.axis("square")
    ax4.axis("off")
    plt.subplots_adjust(left=0.08,right=0.99,top=0.95,bottom=0.15,wspace=0.1,hspace=0.25)
    plt.savefig("figures/objective_layouts_revision.pdf",transparent=True)

    plt.show()



