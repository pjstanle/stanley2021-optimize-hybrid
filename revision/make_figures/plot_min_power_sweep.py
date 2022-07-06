import numpy as np
import matplotlib.pyplot as plt


coe = np.array([51.50323391, 57.23860204, 63.88388862, 70.319386  , 78.6577618 ])
wind = np.array([245000., 245000., 224000., 217000., 196000.])
solar = np.array([257401.04166667, 255256.03298611, 349636.41493056, 379666.53645834,
       538397.17881945])
battery = np.array([125122.18963832, 173998.04496579, 218963.83186706, 265884.65298143,
       312805.4740958 ])

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
plt.gca().set_yticks((50,60,70,80))
plt.gca().set_ylim(48,82)
ax = plt.gca()


plt.subplot(132)
plt.bar(outage,wind/1000,width=1.5,bottom=solar/1000,label="wind")
plt.bar(outage,solar/1000,width=1.5,bottom=0,label="solar")

plt.xlabel("minimum power (MW)", fontsize=8,labelpad=0)
plt.ylabel("capacity (MW)", fontsize=8,labelpad=0)
plt.legend(loc=4,fontsize=8)
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

plt.savefig("figures/sweep_min_power_3.pdf",transparent=True)
plt.show()