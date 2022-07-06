
import numpy as np
import matplotlib.pyplot as plt

from analyze_wind import init_wind_plant
from math import radians



if __name__=="__main__":


    filename = "low_2_43r_116d_88h"
    powercurve_filename = '/Users/astanley/Projects/finished_projects/stanley2021-optimize-hybrid/revision/turbine_data/%s.txt'%filename
    rotor_diameter = 116.0
    hub_height = 88.0
        
    plant = init_wind_plant(hub_height,rotor_diameter,powercurve_filename)

    data = plant.site.wind_resource.data["data"]

    nhours = np.shape(data)[0]
    m0 = np.zeros(nhours)
    m1 = np.zeros(nhours)
    speed = np.zeros(nhours)
    direction = np.zeros(nhours)
    for i in range(nhours):
        m0[i] = data[i][0]
        m1[i] = data[i][1]
        speed[i] = data[i][2]
        direction[i] = data[i][3]

    time = np.arange(nhours)


    ndirs = 36
    dirs = np.linspace(0,360-360/ndirs,ndirs)
    speed_avg = np.zeros(ndirs)
    count = np.zeros(ndirs)

    for i in range(nhours):
        for j in range(ndirs-1):
            if direction[i]%360 < dirs[j+1] and direction[i]%360 > dirs[j]:
                count[j] += 1
                speed_avg[j] += speed[i]
        if direction[i]%360 > dirs[-1]:
            count[-1] += 1
            speed_avg[-1] += speed[i]

    freq = count/nhours
    speed_avg = speed_avg/count

    
    bottom = 0
    width = (2*np.pi) / ndirs

    wd = dirs
    wf = freq
    ws = speed_avg

    # wd -= wd[np.argmax(wf)]
    wd += 270.
    wd +=180./float(ndirs)
    for i in range(ndirs):
        wd[i] = -radians(wd[i])

    nspeeds = 22
    speed_lim = np.linspace(0,21,nspeeds)-0.5
    count = np.zeros(nspeeds)

    for i in range(nhours):
        low = False
        for j in range(nspeeds-1):
            if speed[i] < speed_lim[j+1] and speed[i] >= speed_lim[j]:
                count[j] += 1
                low = True
        if low == False:
            count[-1] += 1
    
    print(max(speed))
    # plt.bar(speed_lim+0.5,count/nhours)
    integral = np.trapz(count, x=speed_lim+0.5)

    print(count[:10])
    plt.figure(figsize=(3.5,2.5))
    plt.fill_between(speed_lim[:11]+0.5,count[:11]/integral,color="C0")
    plt.fill_between(speed_lim[10:]+0.5,count[10:]/integral,color="C1")
    plt.xlabel("wind speed (m/s)",fontsize=8)
    plt.ylabel("probability of occurance",fontsize=8)
    plt.xticks(fontsize=8)
    plt.yticks(fontsize=8)
    plt.gca().set_yticks((0,0.05,0.1))

    print(speed_lim[:11])
    print(speed_lim[10:])

    print(np.trapz(count[10:]/integral,x=speed_lim[10:]+0.5))
    print(np.trapz(count[:11]/integral,x=speed_lim[:11]+0.5))

    plt.text(6,0.04,"79.9%",fontsize=8,horizontalalignment="center")
    plt.text(12,0.01,"21.1%",fontsize=8,horizontalalignment="center")

    plt.tight_layout()

    powers = (speed_lim+0.5)**3
    for i in range(len(powers)):
        if powers[i] > 10**3:
            powers[i] = 10**3

    real_power = sum(powers*count)
    ideal_power = nhours*10**3

    print(real_power)
    print(ideal_power)
    print(real_power/ideal_power)  

    # """wind rose"""
    # plt.figure(figsize=(5,2.5))
    # max_height = max(wf)
    # ax = plt.subplot(121, polar=True)
    # bars = ax.bar(wd, wf*100, width=width, bottom=bottom, color='C0')
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=8)
    # ax.set_yticks((1.5,3,4.5))
    # ax.set_yticklabels(("1.5%","3%","4.5%"),fontsize=8,horizontalalignment="center")
    # ax.set_title("wind frequencies",fontsize=8)
    # """wind speeds"""
    # max_height = max(ws)
    # ax = plt.subplot(122, polar=True)
    # bars = ax.bar(wd, ws, width=width, bottom=bottom, color='C1')
    # ax.set_xticklabels(['E', 'NE', 'N', 'NW', 'W', 'SW', 'S', 'SE'],fontsize=8)
    # ax.set_yticks((5,7,9))
    # ax.set_yticklabels(("5 m/s","7 m/s","9 m/s"),fontsize=8,horizontalalignment="center")
    # ax.set_ylim(0,9)
    # ax.set_title("wind speeds",fontsize=8)

    # plt.subplots_adjust(top=0.78,bottom=0.12,left=0.05,right=0.95,wspace=0.25)



    # plt.figure(figsize=(6,2))
    # ax = plt.subplot(131)
    # ax.plot(time[0:72],direction[0:72])
    # ax.set_xticks((0,12,24,36,48,60,72))
    # plt.xticks(fontsize=8)
    # plt.yticks(fontsize=8)
    # ax.set_xlabel("time (hr)",fontsize=8)
    # ax.set_ylabel("wind direction (degrees)",fontsize=8)

    # ax = plt.subplot(132)
    # ax.plot(time[0:72],speed[0:72],color="C1")
    # ax.set_xticks((0,12,24,36,48,60,72))
    # plt.xticks(fontsize=8)
    # plt.yticks(fontsize=8)
    # ax.set_xlabel("time (hr)",fontsize=8)
    # ax.set_ylabel("wind speed (m/s)",fontsize=8)


    # np.random.seed(seed=1)
    # solar_array = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.05,0.15,0.3,0.5,0.7,0.85,0.95,1.00,1.00,0.95,0.85,0.7,0.5,0.3,0.15,0.05,0.0,0.0])
    # sr_array = np.array([])
    # for i in range(365):
    #     sr_array = np.append(sr_array,solar_array)

    # sr_array = sr_array - np.random.rand(len(sr_array))*0.1
    # for i in range(len(sr_array)):
    #     if sr_array[i] < 0.0:
    #         sr_array[i] = 0.0

    # ax = plt.subplot(133)
    # ax.plot(time[0:72],sr_array[0:72],color="C2")
    # ax.set_xticks((0,12,24,36,48,60,72))
    # plt.xticks(fontsize=8)
    # plt.yticks(fontsize=8)
    # ax.set_xlabel("time (hr)",fontsize=8)
    # ax.set_ylabel("solar capacity factor",fontsize=8)

    # plt.subplots_adjust(top=0.98,bottom=0.2,left=0.1,right=0.99,wspace=0.5)

    # plt.savefig("figures/blue_creek_windrose.pdf",transparent=True)
    # plt.savefig("figures/hourly_resource.pdf",transparent=True)
    # plt.savefig("figures/speeds_pdf.pdf",transparent=True)


    # plt.show()