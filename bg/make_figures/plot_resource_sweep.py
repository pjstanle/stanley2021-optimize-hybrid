import numpy as np
import matplotlib.pyplot as plt

solar = np.array([366.7088572349874,340.4491278426075,304.429782514541,236.16371007188627,107.31302788405742,118.08185503594314])
wind = np.array([28.0,105.0,231.0,266.0,301.0,308.0])
battery = np.array([428.57142857142856,361.9047619047619,361.9047619047619,361.9047619047619,361.9047619047619,361.9047619047619])
coe = np.array([71.99168471923504,65.20042289271726,58.69435155449355,53.0882913748826,48.498805486310786,45.89025527443974])

wsm = np.array([0.6,0.8,1.0,1.2,1.4])

plt.figure(figsize=(6,2.2))

end = 5
plt.subplot(131)
plt.plot(wsm[0:end],coe[0:end],"-o",color="C2")
plt.xlabel("wind speed multiplier", fontsize=8,labelpad=0)
plt.ylabel("COE ($/MWh)", fontsize=8,labelpad=0)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((0.6,0.8,1.0,1.2,1.4))


plt.subplot(132)
plt.bar(wsm[0:end],wind[0:end],width=0.1,bottom=solar[0:end],label="wind")
plt.bar(wsm[0:end],solar[0:end],width=0.1,bottom=0,label="solar")

plt.xlabel("wind speed multiplier", fontsize=8,labelpad=0)
plt.ylabel("capacity (MW)", fontsize=8,labelpad=0)
plt.legend(loc=3,fontsize=8)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((0.6,0.8,1.0,1.2,1.4))

plt.subplot(133)
plt.plot(wsm[0:end],battery[0:end],"-o",color="C3")
plt.xlabel("wind speed multiplier", fontsize=8,labelpad=0)
plt.ylabel("battery capacity (MWh)", fontsize=8,labelpad=0)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.gca().set_xticks((0.6,0.8,1.0,1.2,1.4))


plt.subplots_adjust(bottom=0.2,top=0.99,right=0.99,left=0.08,wspace=0.4)

# plt.savefig("figures/sweep_resource.pdf",transparent=True)
plt.show()