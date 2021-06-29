import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from math import sin, cos, radians
from shapely.geometry import Polygon

if __name__=='__main__':

    
    fig = plt.figure(figsize=[2,2])
    boundary_poly = Polygon(([0,0],[1000,-1000],[1200,1000],[500,1100],[-200,1000])).buffer(2000)

    bx,by = boundary_poly.exterior.coords.xy
    plt.plot(bx,by)
    plt.axis("equal")
    plt.axis("off")

    L = 100
    
    plt.plot([min(bx)-5*L,min(bx)-5*L],[min(by),max(by)],"-k")
    plt.plot([min(bx)-6*L,min(bx)-4*L],[min(by),min(by)],"-k")
    plt.plot([min(bx)-6*L,min(bx)-4*L],[max(by),max(by)],"-k")

    plt.plot([min(bx),max(bx)],[min(by)-5*L,min(by)-5*L],"-k")
    plt.plot([min(bx),min(bx)],[min(by)-6*L,min(by)-4*L],"-k")
    plt.plot([max(bx),max(bx)],[min(by)-6*L,min(by)-4*L],"-k")

    height = 6.1
    width = 5.4
    area = 26.1

    plt.text(min(bx)+(width*1000)/2.0,min(by)-7*L,"%s km"%width,horizontalalignment="center",verticalalignment="top",fontsize=8)
    plt.text(min(bx)-7*L,min(by)+(height*1000)/2.0,"%s km"%height,horizontalalignment="right",verticalalignment="center",fontsize=8,rotation=90)   
    plt.text(min(bx)+(width*1000)/2.0,min(by)+(height*1000)/2.0,r"%s km$^2$"%area,horizontalalignment="center",verticalalignment="center",fontsize=8)
    


    plt.subplots_adjust(top = 0.99, bottom = 0.1, right = 0.99, left = 0.1)
    plt.savefig("figures/boundary.pdf",transparent=True)
    plt.show()
