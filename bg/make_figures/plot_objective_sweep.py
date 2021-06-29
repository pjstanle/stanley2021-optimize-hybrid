import numpy as np
import matplotlib.pyplot as plt







# pro min
# solar capacity:  309.96486946935073
# wind capacity:  231.0
# battery capacity:  361.9047619047619
# aep:  -1564450.5580223498
# coe:  58.74858043564813
# pro:  -17.60228961599387


# coe no min
# solar capacity:  303.89857674035636
# wind capacity:  0.0
# battery capacity:  0.0
# aep:  -911897.8293194273
# coe:  34.961294849676946
# pro:  -31.951719168743036

plt.figure(figsize=(4.5,4))

plt.subplot(221)
plt.bar([1],[303.89857674035636/562.7251013488436],width=0.5,color="C1")
plt.text(1,303.89857674035636/2/562.7251013488436,"304 MW",fontsize=6,horizontalalignment="center",verticalalignment="center")
plt.bar([1],[0/562.7251013488436],width=0.5,color="C0",bottom=303.89857674035636/562.7251013488436)
plt.plot([2],[0.0/361.9047619047619],"o",color="C3")
plt.text(2,0.0/361.9047619047619+0.1,"0 MWh",fontsize=6,horizontalalignment="center",verticalalignment="center")
plt.plot([3],[34.961294849676946/58.74858043564813],"o",color="C2")
plt.text(3,34.961294849676946/58.74858043564813+0.1,"34.96 $/MWh",fontsize=6,horizontalalignment="center",verticalalignment="center")
plt.plot([4],[31.951719168743036/47.364312588801475],"o",color="black")
plt.text(4,31.951719168743036/47.364312588801475+0.1,"32.0 $MM",fontsize=6,horizontalalignment="center",verticalalignment="center")

plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)

plt.gca().set_yticks(())
plt.ylim(0,1.05)
plt.gca().set_xticks((1,2,3,4))
plt.gca().set_xticklabels(("generation","battery","COE","profit"),fontsize=8,rotation=30)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title("minimize COE",fontsize=8)
plt.ylabel("no minimum\npower constraint",fontsize=8,rotation=0,labelpad=40)


# coe min
# solar capacity:  309.96486946935073
# wind capacity:  231.0
# battery capacity:  361.9047619047619
# aep:  -1564450.5580223498
# coe:  58.74858043564813
# pro:  -17.60228961599387
plt.subplot(223)
plt.bar([1],[309.96486946935073/562.7251013488436],width=0.5,color="C1")
plt.text(1,309.96486946935073/2/562.7251013488436,"310 MW",fontsize=6,horizontalalignment="center",verticalalignment="center")
plt.bar([1],[231.0/562.7251013488436],width=0.5,color="C0",bottom=309.96486946935073/562.7251013488436)
plt.text(1,(309.96486946935073+231/2)/562.7251013488436,"231 MW",fontsize=6,horizontalalignment="center",verticalalignment="center")
plt.plot([2],[361.9047619047619/361.9047619047619],"o",color="C3")
plt.plot([3],[58.74858043564813/58.74858043564813],"o",color="C2")
plt.plot([4],[17.60228961599387/47.364312588801475],"o",color="black")
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.gca().set_yticks(())
plt.ylim(0,1.05)
plt.gca().set_xticks((1,2,3,4))
plt.gca().set_xticklabels(("generation","battery","COE","profit"),fontsize=8,rotation=30)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.ylabel("30 MW minimum\npower constraint",fontsize=8,rotation=0,labelpad=40)


# pro no min
# solar capacity:  324.7251013488436
# wind capacity:  238.0
# battery capacity:  0.0
# aep:  -1599856.7645964464
# coe:  40.39465429847475
# pro:  -47.364312588801475

plt.subplot(222)
plt.bar([1],[324.7251013488436/562.7251013488436],width=0.5,color="C1",label="solar")
plt.text(1,324.7251013488436/2/562.7251013488436,"325 MW",fontsize=6,horizontalalignment="center",verticalalignment="center")
plt.bar([1],[238.0/562.7251013488436],width=0.5,color="C0",bottom=324.7251013488436/562.7251013488436,label="wind")
plt.text(1,(324.7251013488436+238.0/2)/562.7251013488436,"238 MW",fontsize=6,horizontalalignment="center",verticalalignment="center")
plt.plot([2],[0.0/361.9047619047619],"o",color="C3")
plt.plot([3],[40.39465429847475/58.74858043564813],"o",color="C2")
plt.plot([4],[47.364312588801475/47.364312588801475],"o",color="black")
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.gca().set_yticks(())
plt.ylim(0,1.05)
plt.gca().set_xticks((1,2,3,4))
plt.gca().set_xticklabels(("generation","battery","COE","profit"),fontsize=8,rotation=30)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)
plt.title("maximize profit",fontsize=8)
plt.legend(fontsize=8,loc=4)


plt.subplot(224)
plt.bar([1],[309.96486946935073/562.7251013488436],width=0.5,color="C1")
plt.text(1,309.96486946935073/2/562.7251013488436,"310 MW",fontsize=6,horizontalalignment="center",verticalalignment="center")
plt.bar([1],[231.0/562.7251013488436],width=0.5,color="C0",bottom=309.96486946935073/562.7251013488436)
plt.text(1,(309.96486946935073+231.0/2)/562.7251013488436,"231 MW",fontsize=6,horizontalalignment="center",verticalalignment="center")
plt.plot([2],[361.9047619047619/361.9047619047619],"o",color="C3")
plt.plot([3],[58.74858043564813/58.74858043564813],"o",color="C2")
plt.plot([4],[17.60228961599387/47.364312588801475],"o",color="black")
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.gca().set_yticks(())
plt.ylim(0,1.05)
plt.gca().set_xticks((1,2,3,4))
plt.gca().set_xticklabels(("generation","battery","COE","profit"),fontsize=8,rotation=30)
plt.xticks(fontsize=8)
plt.yticks(fontsize=8)


plt.subplots_adjust(left=0.25,bottom=0.15,top=0.9,hspace=0.6,right=0.97)

plt.show()