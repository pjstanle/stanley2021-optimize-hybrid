import numpy as np
import matplotlib.pyplot as plt
from numpy.core.numeric import full
from plotting_functions import plot_poly, plot_turbines
from shapely.geometry import Polygon


if __name__=="__main__":

    corner = 1.0
    boundary_polygon = Polygon(((-corner,-corner),(-corner,corner),(corner,corner),(corner,-corner)))

    plt.figure(figsize=(6,2))
    ax1 = plt.subplot(141)
    ax2 = plt.subplot(142)
    ax3 = plt.subplot(143)
    ax4 = plt.subplot(144)

    plot_poly(boundary_polygon,ax=ax1,color="black",fill=False)
    plot_poly(boundary_polygon,ax=ax2,color="black",fill=False)
    plot_poly(boundary_polygon,ax=ax3,color="black",fill=False)
    plot_poly(boundary_polygon,ax=ax4,color="black",fill=False)


    ax1.axis("square")
    ax2.axis("square")
    ax3.axis("square")
    ax4.axis("square")

    ax1.axis("off")
    ax2.axis("off")
    ax3.axis("off")
    ax4.axis("off")

    np.random.seed(5)
    N = 6
    p = np.random.rand(2*N)*2*corner-corner
    px = p[0:N]
    py = p[N:2*N]

    ax1.scatter(px,py,color="C0")
    ax1.set_title("parent population",fontsize=8)



    pairs = np.arange(N)
    np.random.seed(2)
    np.random.shuffle(pairs)
    np.random.seed(6)

    full_x = np.copy(px)
    full_y = np.copy(py)
    for i in range(int(N/2)):
        x = [px[pairs[2*i]],px[pairs[2*i+1]]]
        y = [py[pairs[2*i]],py[pairs[2*i+1]]]
        ax2.plot(x,y,color="C%s"%i)
        ax2.scatter(x,y,color="C%s"%i)

        dx = abs(px[pairs[2*i]] - px[pairs[2*i+1]])
        new_x = (px[pairs[2*i]] + px[pairs[2*i+1]])/2 + dx*(np.random.rand(2)-0.5) 
        dy = abs(py[pairs[2*i]] - py[pairs[2*i+1]])
        new_y = (py[pairs[2*i]] + py[pairs[2*i+1]])/2 + dy*(np.random.rand(2)-0.5) 
        ax2.scatter(new_x,new_y,color="C%s"%i,facecolors='none')

        full_x = np.append(full_x,new_x)
        full_y = np.append(full_y,new_y)

    ax2.set_title("crossover",fontsize=8)

    d = 0.1
    for i in range(N*2):
        x = full_x[i] + d*(np.random.rand()*2-1)
        y = full_y[i] + d*(np.random.rand()*2-1)

        ax3.scatter(full_x[i],full_y[i],color="C0",facecolors='none')
        ax3.plot([full_x[i],x],[full_y[i],y],"--",color="black",linewidth=0.5)
        ax3.scatter(x,y,color="C0")

        full_x[i] = x
        full_y[i] = y

    ax3.set_title("mutation", fontsize=8)


    indices = [0,2,4,6,8,10]
    ax4.scatter(full_x[indices],full_y[indices],s=100,color="black",facecolors='none')
    ax4.scatter(full_x,full_y,color="C0")
    ax4.set_title("selection",fontsize=8)
    
    plt.subplots_adjust(left=0.0,right=1.0,top=0.95,bottom=0.01,wspace=0.05)
    plt.savefig("figures/ga_overview.pdf",transparent=True)
    plt.show()



