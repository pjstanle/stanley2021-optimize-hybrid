import numpy as np
import matplotlib.pyplot as plt


coe = np.array([214.50653627, 100.60608828,  63.44846289,  49.13917911,
        41.32051758])
wind = np.array([147000., 196000., 238000., 245000., 245000.])
solar = np.array([143715.58159722, 538397.17881945, 285286.15451389, 293866.18923611,
       240240.97222222])
battery = np.array([322580.64516129, 250244.37927664, 226783.96871945, 215053.76344086,
       187683.28445748])

# coe = np.array([59.20986284, 49.8182307 , 34.54051256, 27.8541512 , 24.63042146])
# wind = np.array([     0.,  98000., 112000., 112000., 140000.])
# solar = np.array([516947.09201389,      0.        ,      0.        ,      0.        ,
#             0.        ])
# battery = np.array([0., 0., 0., 0., 0.])

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