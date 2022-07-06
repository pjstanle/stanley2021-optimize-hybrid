import numpy as np
import matplotlib.pyplot as plt


solar = np.array([326.11583366, 325.21031087, 332.10521729, 325.70441194,
       309.96486947, 309.96486947,408.0680985447498,335.7952752584633,345.0204201831463,381.92099988187863])
wind = np.array([196., 196., 196., 196., 224., 252.,231.0,252.0,252.0,224.0])
battery = np.array([304.76190476, 304.76190476, 304.76190476, 304.76190476,
       361.9047619 , 457.14285714,571.4285714285714,942.8571428571429,1085.7142857142858,1657.142857142857])
coe = np.array([55.80524101, 55.98073369, 56.03972098, 55.98078909, 58.76974033,
       63.49768293,71.63733544389433,87.9565189249901,95.09106717792392,125.88722848676245])

outage = np.array([0,3,6,9,12,15,18,24,36,48])

# plt.figure(figsize=(6,2.2))
plt.figure(figsize=(2.2,5.5))

end = 8
plt.subplot(311)
plt.plot(outage[0:end],coe[0:end],"-o",color="C2")
# plt.xlabel("outage duration (hr)", fontsize=12,labelpad=0)
plt.ylabel("COE ($/MWh)", fontsize=12,labelpad=0)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.gca().set_xticks((0,6,12,18,24))


plt.subplot(312)
plt.bar(outage[0:end],wind[0:end],width=1.5,bottom=solar[0:end],label="wind")
plt.bar(outage[0:end],solar[0:end],width=1.5,bottom=0,label="solar")

# plt.xlabel("outage duration (hr)", fontsize=12,labelpad=0)
plt.ylabel("capacity (MW)", fontsize=12,labelpad=0)
plt.legend(loc=3,fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.gca().set_xticks((0,6,12,18,24))

plt.subplot(313)
plt.plot(outage[0:end],battery[0:end],"-o",color="C3")
plt.xlabel("outage duration (hr)", fontsize=12,labelpad=0)
plt.ylabel("battery capacity (MWh)", fontsize=12,labelpad=0)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.gca().set_xticks((0,6,12,18,24))


plt.subplots_adjust(bottom=0.1,top=0.99,right=0.88,left=0.28,wspace=0.4)

# plt.savefig("figures/sweep_outage.pdf",transparent=True)
# plt.savefig("figures/sweep_outage_ab.pdf",transparent=True)
plt.show()