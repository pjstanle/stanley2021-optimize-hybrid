import numpy as np
import matplotlib.pyplot as plt


coe = np.array([61.46484866, 62.74827445, 62.78957669, 63.60480828, 65.64673779])
wind = np.array([245000., 245000., 238000., 231000., 224000.])
solar = np.array([298156.20659722, 255256.03298611, 311026.25868056, 326041.31944445,
       296011.19791667])
battery = np.array([222873.90029326, 220918.86608016, 219941.34897361, 218963.83186706,
       222873.90029326])

outage = np.array([6,8,10,12,14])

plt.figure(figsize=(6,1.8))

from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

plt.subplot(131)
plt.plot(outage,coe,"-o",color="C2")
plt.xlabel("minimum power (MW)", fontsize=8,labelpad=0)
plt.ylabel("COE ($/MWh)", fontsize=8,labelpad=0)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((6,8,10,12,14))
ax = plt.gca()


plt.subplot(132)
plt.bar(outage,wind/1000,width=1.5,bottom=solar/1000,label="wind")
plt.bar(outage,solar/1000,width=1.5,bottom=0,label="solar")

plt.xlabel("minimum power (MW)", fontsize=8,labelpad=0)
plt.ylabel("capacity (MW)", fontsize=8,labelpad=0)
plt.legend(loc=3,fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((6,8,10,12,14))

plt.subplot(133)
plt.plot(outage,battery/1000,"-o",color="C3")
plt.xlabel("minimum power (MW)", fontsize=8,labelpad=0)
plt.ylabel("battery capacity (MWh)", fontsize=8,labelpad=0)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((6,8,10,12,14))


plt.subplots_adjust(bottom=0.2,top=0.99,right=0.99,left=0.08,wspace=0.6)

# plt.savefig("figures/sweep_min_power_2.pdf",transparent=True)
plt.show()