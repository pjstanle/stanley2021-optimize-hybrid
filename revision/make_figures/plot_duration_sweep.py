import numpy as np
import matplotlib.pyplot as plt


# coe = np.array([64.31458222, 63.28776946, 63.92008526, 63.8176058 , 64.78892647,
#        65.73606573, 68.676496  , 70.13886517, 75.22409722, 75.24213535,
#        78.02908782, 80.78332876, 83.89300411])
# wind = np.array([224000., 231000., 245000., 231000., 266000., 294000., 287000.,
#        301000., 301000., 308000., 308000., 301000., 301000.])
# solar = np.array([313171.26736111, 328186.328125  , 296011.19791667, 300301.21527778,
#        255256.03298611, 169455.68576389, 120120.48611111,  83655.33854167,
#         42900.17361111,  25740.10416667,  25740.10416667,  53625.21701389,
#         42900.17361111])
# battery = np.array([222873.90029326, 218963.83186706, 222873.90029326, 221896.38318671,
#        250244.37927664, 283479.96089932, 304985.3372434 , 330400.78201369,
#        391006.84261975, 391006.84261975, 420332.35581623, 450635.38611926,
#        480938.41642229])

coe = np.array([64.31458222, 63.28776946, 63.92008526, 63.8176058 , 64.78892647,
       65.73606573, 68.676496  , 70.13886517, 75.24213535,
       78.02908782, 80.78332876, 83.89300411])
wind = np.array([224000., 231000., 245000., 231000., 266000., 294000., 287000.,
       301000., 308000., 308000., 301000., 301000.])
solar = np.array([313171.26736111, 328186.328125  , 296011.19791667, 300301.21527778,
       255256.03298611, 169455.68576389, 120120.48611111,  83655.33854167,
       25740.10416667,  25740.10416667,  53625.21701389,
        42900.17361111])
battery = np.array([222873.90029326, 218963.83186706, 222873.90029326, 221896.38318671,
       250244.37927664, 283479.96089932, 304985.3372434 , 330400.78201369,
       391006.84261975, 420332.35581623, 450635.38611926,
       480938.41642229])

outage = np.array([0,6,12,18,24,27,30,33,39,42,45,48])

plt.figure(figsize=(6,1.8))

from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

plt.subplot(131)
plt.plot(outage,coe,"-o",color="C2")
plt.xlabel("outage duration (hr)", fontsize=8,labelpad=0)
plt.ylabel("COE ($/MWh)", fontsize=8,labelpad=0)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((0,8,16,24,32,40,48))
ax = plt.gca()

# axins = ax.inset_axes([0.25,0.5,0.45,0.4])
# axins.plot(outage[0:-3],coe[0:-3],"-o",markersize=4,color="C2")
# axins.set_xticks((8,12,16,20,24))
# axins.set_xticklabels(("8","12","16","20","24"),fontsize=6)

# axins.set_yticks((100.4,100.5,100.6))
# axins.set_yticklabels(("100.4","100.5","100.6"),fontsize=6)

# mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.75")


plt.subplot(132)
plt.bar(outage,wind/1000,width=1.5,bottom=solar/1000,label="wind")
plt.bar(outage,solar/1000,width=1.5,bottom=0,label="solar")

plt.xlabel("outage duration (hr)", fontsize=8,labelpad=0)
plt.ylabel("capacity (MW)", fontsize=8,labelpad=0)
plt.legend(loc=3,fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((0,8,16,24,32,40,48))

plt.subplot(133)
plt.plot(outage,battery/1000,"-o",color="C3")
plt.xlabel("outage duration (hr)", fontsize=8,labelpad=0)
plt.ylabel("battery capacity (MWh)", fontsize=8,labelpad=0)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((0,8,16,24,32,40,48))


plt.subplots_adjust(bottom=0.2,top=0.99,right=0.99,left=0.08,wspace=0.6)

plt.savefig("figures/sweep_outage_3.pdf",transparent=True)
plt.show()