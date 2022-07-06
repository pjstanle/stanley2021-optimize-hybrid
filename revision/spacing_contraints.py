import numpy as np

class WindSpacingConstraints:
    
    
    def __init__(self, turbine_x, turbine_y):
        # wind
        self.turbine_x = turbine_x
        self.turbine_y = turbine_y


    def get_wind_wind_distance(self):

        #calculate the spacing between each turbine and every other turbine (without repeating)
        nturbs = len(self.turbine_x)

        if nturbs == 1:
            spacing = np.array([1E6])
        
        else:
            npairs = int((nturbs*(nturbs-1))/2)
            spacing = np.zeros(npairs)

            ind = 0
            for i in range(nturbs):
                for j in range(i,nturbs):
                    if i != j:
                        spacing[ind] = np.sqrt((self.turbine_x[i]-self.turbine_x[j])**2+(self.turbine_y[i]-self.turbine_y[j])**2)
                        ind += 1

        return spacing

        
if __name__=="__main__":

    N = 100
    turbine_x = np.arange(N)
    turbine_y = np.zeros(N)

    spacing = WindSpacingConstraints(turbine_x, turbine_y)
    s = spacing.get_wind_wind_distance()
    print(s)