# implement a simple dispatch model

import numpy as np

class SimpleDispatch():

    def __init__(self, Nt, min_power):

        # length of simulation
        self.Nt = Nt
        self.min_power = min_power

        # amount of curtailment experienced by plant
        self.plant_power = np.zeros(self.Nt)
        
        # size of battery (MWh)
        self.battery_size = 0
        

    def run(self):

        battery_SOC = np.zeros(self.Nt)
        battery_discharged = np.zeros(self.Nt)
        battery_charged = np.zeros(self.Nt)
        battery_SOC[0] = self.battery_size
        for i in range(self.Nt):

            # should you charge
            if self.plant_power[i] > self.min_power:
                excess = self.plant_power[i] - self.min_power
                # if i == 0:
                #     battery_SOC[i] = np.min([excess, self.battery_size])
                #     battery_charged[i] = battery_SOC[i]
                if i > 0:
                    if battery_SOC[i-1] < self.battery_size:
                        add_gen = np.min([excess, self.battery_size])
                        battery_SOC[i] = np.min([battery_SOC[i-1] + add_gen, self.battery_size])
                        battery_charged[i] = battery_SOC[i] - battery_SOC[i-1]
                    else:
                        battery_SOC[i] = battery_SOC[i - 1]

            # should you discharge
            else:
                power_needed = self.min_power - self.plant_power[i]
                if i > 0:
                    if battery_SOC[i-1] > 0:
                        battery_discharged[i] = np.min([power_needed, battery_SOC[i-1], self.battery_size])
                        battery_SOC[i] = battery_SOC[i-1] - battery_discharged[i]

        return battery_charged, battery_discharged, battery_SOC



