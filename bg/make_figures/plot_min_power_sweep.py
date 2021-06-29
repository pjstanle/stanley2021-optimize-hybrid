import numpy as np
import matplotlib.pyplot as plt

solar = np.array([303.78402446, 299.81721005, 297.04966657, 309.96486947,
       324.72510135, 324.72510135])
wind = np.array([  0., 119., 196., 224., 231., 238.])
battery = np.array([  0.        , 126.98412698, 253.96825397, 361.9047619 ,
       507.93650794, 603.17460317])
coe = np.array([34.93517855, 45.50729925, 52.96760687, 58.76974033, 66.24130686,
       70.96408881])

min_power = np.array([0,10,20,30,40,50])

plt.figure(figsize=(6,2.2))

plt.subplot(131)
plt.plot(min_power,coe,"-o",color="C2")
plt.xlabel("minimum power (MW)", fontsize=8,labelpad=0)
plt.ylabel("COE ($/MWh)", fontsize=8,labelpad=0)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((0,20,40))



plt.subplot(132)
plt.bar(min_power,wind,width=5,bottom=solar,label="wind")
plt.bar(min_power,solar,width=5,bottom=0,label="solar")

plt.xlabel("minimum power (MW)", fontsize=8,labelpad=0)
plt.ylabel("capacity (MW)", fontsize=8,labelpad=0)
plt.legend(loc=3,fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((0,20,40))

plt.subplot(133)
plt.plot(min_power,battery,"-o",color="C3")
plt.xlabel("minimum power (MW)", fontsize=8,labelpad=0)
plt.ylabel("battery capacity (MWh)", fontsize=8,labelpad=0)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((0,20,40))


plt.subplots_adjust(bottom=0.2,top=0.99,right=0.99,left=0.08,wspace=0.4)

plt.savefig("figures/sweep_min_power.pdf",transparent=True)
plt.show()