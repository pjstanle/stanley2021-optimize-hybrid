import numpy as np
import matplotlib.pyplot as plt

solar = np.array([332.10521729, 318.38432877, 294.2821231 , 295.20463759,
       249.07891297, 236.16371007, 125.46197098])
wind = np.array([175., 189., 224., 231., 273., 259., 287.])
battery = np.array([361.9047619, 361.9047619, 361.9047619, 361.9047619, 361.9047619,
       361.9047619, 361.9047619])
coe = np.array([52.09821553, 56.33124665, 59.61596498, 63.06540833, 65.94567549,
       69.07673862, 70.98603776])

solar_multiplier = np.array([0.4,0.6,0.8,1.0,1.2,1.4,1.6])

plt.figure(figsize=(6,2.2))

plt.subplot(131)
plt.plot(solar_multiplier,coe,"-o",color="C2")
plt.xlabel("solar cost multiplier", fontsize=8,labelpad=0)
plt.ylabel("COE ($/MWh)", fontsize=8,labelpad=0)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((0.4,0.8,1.2,1.6))


plt.subplot(132)
plt.bar(solar_multiplier,wind,width=0.1,bottom=solar,label="wind")
plt.bar(solar_multiplier,solar,width=0.1,bottom=0,label="solar")

plt.xlabel("solar cost multiplier", fontsize=8,labelpad=0)
plt.ylabel("capacity (MW)", fontsize=8,labelpad=0)
plt.legend(loc=3,fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((0.4,0.8,1.2,1.6))

plt.subplot(133)
plt.plot(solar_multiplier,battery,"-o",color="C3")
plt.xlabel("solar cost multiplier", fontsize=8,labelpad=0)
plt.ylabel("battery capacity (MWh)", fontsize=8,labelpad=0)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((0.4,0.8,1.2,1.6))


plt.subplots_adjust(bottom=0.2,top=0.99,right=0.99,left=0.08,wspace=0.4)

plt.savefig("figures/sweep_cost.pdf",transparent=True)
plt.show()