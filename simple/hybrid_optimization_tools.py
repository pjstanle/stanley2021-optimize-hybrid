import numpy as np
import pyoptsparse

def hybrid_greedy(obj_function,nlocs,ntech,initialize=0): 

    plant_array = np.zeros(nlocs)
    plant_solution = 1E6
    while plant_solution > 1000.0:
        print("searching")
        plant_array[0:initialize] = np.random.randint(1,high=(ntech+1),size=initialize)
        np.random.shuffle(plant_array)
        plant_solution = obj_function(plant_array)

    temp_array = np.zeros_like(plant_array)
    best_array = np.zeros_like(plant_array)
    
    converged = False
    while converged == False:
        best_solution = plant_solution
        best_array[:] = plant_array[:]
        for i in range(nlocs):
            for j in range(ntech+1):
                temp_array[:] = plant_array[:]
                temp_array[i] = j
                temp_solution = obj_function(temp_array)
                if temp_solution < best_solution:
                    best_solution = temp_solution
                    best_array[:] = temp_array[:]
        
        if best_solution == plant_solution:
            converged = True
            solution = plant_solution
            array = plant_array
        else:
            plant_solution = best_solution
            plant_array[:] = best_array[:]
            print(plant_solution)
    
    return solution, array


def hybrid_sweep(obj_function,nrows,ncols,ntech): 

    nlocs = nrows*ncols
    plant_array = np.zeros(nlocs)

    # plant_array = np.array([0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
    #    1.000, 0.000, 1.000, 0.000, 0.000, 2.000, 0.000, 0.000, 0.000,
    #    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.000,
    #    0.000, 0.000, 1.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.000,
    #    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.000,
    #    0.000, 0.000, 2.000, 0.000, 1.000, 0.000, 0.000, 0.000, 0.000,
    #    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.000, 0.000,
    #    1.000, 0.000, 0.000, 0.000, 0.000, 2.000, 0.000, 0.000, 0.000,
    #    1.000, 0.000, 0.000, 2.000, 2.000, 0.000, 2.000, 2.000, 0.000,
    #    0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
    #    0.000, 0.000, 0.000, 2.000, 0.000, 1.000, 0.000, 0.000, 2.000,
    #    0.000])
    plant_solution = obj_function(plant_array)
    
    # TODO randomly initialize starting array
    # while plant_solution > 1000.0:
    #     print("searching")
    #     plant_array[0:initialize] = np.random.randint(1,high=(ntech+1),size=initialize)
    #     np.random.shuffle(plant_array)
    #     plant_solution = obj_function(plant_array)

    converged = False
    phase = np.array(["search","switch_x","switch_y"])
    phase_iter = 1 # start with the search phase
    
    converged_counter = 0

    while converged == False:
        phase_iter = (phase_iter+1)%len(phase)
        current_phase = phase[phase_iter]
        # print("plant array: ", plant_array)
        # print("plant solution: ", plant_solution)

        # this is the search phase
        # sweep through every point, and try every technology at this point
        # if the solution is better, keep the new change
        
        if current_phase == "search":
            print("search")
            order = np.arange(nlocs)
            np.random.shuffle(order)
            for i in range(nlocs):
                temp_array = np.zeros_like(plant_array)
                temp_array[:] = plant_array[:]
                for j in range(ntech):
                    temp_array[order[i]] = (temp_array[order[i]]+1)%(ntech+1)
                    temp_solution = obj_function(temp_array)
                    if temp_solution < plant_solution:
                        plant_solution = temp_solution
                        plant_array[order[i]] = temp_array[order[i]]
                        converged_counter = 0
                        print(plant_solution)
                    else:
                        converged_counter += 1
                
                if converged_counter > (len(phase)+1)*nlocs:
                    converged = True
                    break
                    
        
        # elif current_phase == "switch_x":
        #     # print("switch x")

        #     order = np.arange(nlocs)
        #     np.random.shuffle(order)

        #     for y in range(nrows):
        #         for x in range(ncols-1):
        #             temp_array = np.zeros_like(plant_array)
        #             temp_array[:] = plant_array[:]
        #             index = y*ncols + x
        #             if temp_array[index] != temp_array[index+1]:
        #                 temp_array[index],temp_array[index+1] = temp_array[index+1],temp_array[index]
        #                 temp_solution = obj_function(temp_array)

        #                 if temp_solution < plant_solution:
        #                     plant_solution = temp_solution
        #                     plant_array[i] = temp_array[i]
        #                     converged_counter = 0
        #                 else:
        #                     converged_counter += 1

        #             else:
        #                 converged_counter += 1

        #     if converged_counter > (len(phase)+1)*nlocs:
        #             converged = True
        #             break


        # elif current_phase == "switch_y":
        #     # print("switch y")
        #     for y in range(nrows-1):
        #         for x in range(ncols):
        #             temp_array = np.zeros_like(plant_array)
        #             temp_array[:] = plant_array[:]
        #             index = y*ncols + x
        #             if temp_array[index] != temp_array[index+1]:
        #                 temp_array[index],temp_array[index+1] = temp_array[index+1],temp_array[index]
        #                 temp_solution = obj_function(temp_array)

        #                 if temp_solution < plant_solution:
        #                     plant_solution = temp_solution
        #                     plant_array[i] = temp_array[i]
        #                     converged_counter = 0
        #                 else:
        #                     converged_counter += 1

        #             else:
        #                 converged_counter += 1

        #     if converged_counter > (len(phase)+1)*nlocs:
        #             converged = True
        #             break
        elif current_phase == "switch_x":
            print("switch x")
            order = np.arange(nlocs)
            np.random.shuffle(order)
            for i in range(nlocs):
                temp_array = np.zeros_like(plant_array)
                temp_array[:] = plant_array[:]
                try:
                    if temp_array[order[i]] != temp_array[order[i]+1]:
                        temp_array[order[i]],temp_array[order[i]+1] = temp_array[order[i]+1],temp_array[order[i]]
                        temp_solution = obj_function(temp_array)
                    
                        if temp_solution < plant_solution:
                            plant_solution = temp_solution
                            plant_array[i] = temp_array[i]
                            converged_counter = 0
                            print(plant_solution)
                        else:
                            converged_counter += 1
                    
                    else:
                        converged_counter += 1

                except:
                    converged_counter += 1

                if converged_counter > (len(phase)+1)*nlocs:
                        converged = True
                        break


        elif current_phase == "switch_y":
            print("switch y")
            order = np.arange(nlocs)
            np.random.shuffle(order)
            for i in range(nlocs):
                temp_array = np.zeros_like(plant_array)
                temp_array[:] = plant_array[:]
                try:
                    if temp_array[order[i]] != temp_array[order[i]+ncols]:
                        temp_array[order[i]],temp_array[order[i]+ncols] = temp_array[order[i]+ncols],temp_array[order[i]]
                        temp_solution = obj_function(temp_array)
                    
                        if temp_solution < plant_solution:
                            plant_solution = temp_solution
                            plant_array[i] = temp_array[i]
                            converged_counter = 0
                            print(plant_solution)
                        else:
                            converged_counter += 1
                    
                    else:
                        converged_counter += 1

                except:
                    converged_counter += 1

                if converged_counter > (len(phase)+1)*nlocs:
                        converged = True
                        break

        
    return plant_solution, plant_array


def obj_func_bat(input_dict):

    # calculate the wind farm AEP as a function of the grid design variables
    global input_plant_array
    global full_obj_func

    global batt_scale

    funcs = {}
    fail = False

    # objective
    battery_storage = input_dict["battery_storage"]*batt_scale
    input_array = np.append(input_plant_array,battery_storage)
    obj = full_obj_func(input_array)
    funcs["obj"] = obj/1E3

    return funcs, fail


def opt_battery(obj_function,plant_array,old_battery_storage):
    
    global input_plant_array
    global full_obj_func

    global batt_scale

    batt_scale = 1.2E2

    full_obj_func = obj_function
    input_plant_array = plant_array
    
    optProb = pyoptsparse.Optimization("BatteryOpt",obj_func_bat)
    optProb.addVar("battery_storage",type="c",lower=0.0,upper=None,value=old_battery_storage/batt_scale)

    optProb.addObj("obj")
    optimize = pyoptsparse.SNOPT()
    solution = optimize(optProb,sens="FD")
    
    new_battery_storage = solution.getDVs()["battery_storage"]*batt_scale

    return new_battery_storage


def hybrid_sweep_w_battery(obj_function,nrows,ncols,ntech,start_battery=0.0): 

    global plant_array
    global bat_obj_function

    bat_obj_function = obj_function

    nlocs = nrows*ncols
    plant_array = np.zeros(nlocs)
    battery_storage = start_battery
    input_array = np.append(plant_array,battery_storage)
    plant_solution = obj_function(input_array)
    
    converged = False
    phase = np.array(["search","switch_x","switch_y"])
    phase_iter = 1 # start with the search phase
    
    converged_counter = 0

    while converged == False:
        phase_iter = (phase_iter+1)%len(phase)
        current_phase = phase[phase_iter]

        if current_phase == "search":
            print("search")
            order = np.arange(nlocs)
            np.random.shuffle(order)
            for i in range(nlocs):
                temp_array = np.zeros_like(plant_array)
                temp_array[:] = plant_array[:]
                for j in range(ntech):
                    temp_array[order[i]] = (temp_array[order[i]]+1)%(ntech+1)
                    temp_input = np.append(temp_array,battery_storage)
                    temp_solution = obj_function(temp_input)
                    if temp_solution < plant_solution:
                        plant_solution = temp_solution
                        plant_array[order[i]] = temp_array[order[i]]
                        converged_counter = 0
                        print(plant_solution)
                    else:
                        converged_counter += 1
                
                if converged_counter > (len(phase)+1)*nlocs:
                    converged = True
                    break
                    
        elif current_phase == "switch_x":
            print("switch x")
            order = np.arange(nlocs)
            np.random.shuffle(order)
            for i in range(nlocs):
                temp_array = np.zeros_like(plant_array)
                temp_array[:] = plant_array[:]
                try:
                    if temp_array[order[i]] != temp_array[order[i]+1]:
                        temp_array[order[i]],temp_array[order[i]+1] = temp_array[order[i]+1],temp_array[order[i]]
                        temp_input = np.append(temp_array,battery_storage)
                        temp_solution = obj_function(temp_input)
                    
                        if temp_solution < plant_solution:
                            plant_solution = temp_solution
                            plant_array[i] = temp_array[i]
                            converged_counter = 0
                            print(plant_solution)
                        else:
                            converged_counter += 1
                    
                    else:
                        converged_counter += 1

                except:
                    converged_counter += 1

                if converged_counter > (len(phase)+1)*nlocs:
                        converged = True
                        break

        elif current_phase == "switch_y":
            print("switch y")
            order = np.arange(nlocs)
            np.random.shuffle(order)
            for i in range(nlocs):
                temp_array = np.zeros_like(plant_array)
                temp_array[:] = plant_array[:]
                try:
                    if temp_array[order[i]] != temp_array[order[i]+ncols]:
                        temp_array[order[i]],temp_array[order[i]+ncols] = temp_array[order[i]+ncols],temp_array[order[i]]
                        temp_input = np.append(temp_array,battery_storage)
                        temp_solution = obj_function(temp_input)
                    
                        if temp_solution < plant_solution:
                            plant_solution = temp_solution
                            plant_array[i] = temp_array[i]
                            converged_counter = 0
                            print(plant_solution)
                        else:
                            converged_counter += 1
                    
                    else:
                        converged_counter += 1

                except:
                    converged_counter += 1

                if converged_counter > (len(phase)+1)*nlocs:
                        converged = True
                        break

        if current_phase == "search":
            print("old battery storage: ", battery_storage)
            battery_storage = opt_battery(obj_function,plant_array,battery_storage)
            print("new battery storage: ", battery_storage)

    battery_storage = opt_battery(obj_function,plant_array,battery_storage)
    temp_input = np.append(plant_array,battery_storage)
    plant_solution = obj_function(temp_input)

    return plant_solution, plant_array, battery_storage


if __name__=="__main__":
    def obj_function(x):
        sol = 0.0
        for i in range(len(x)):
            sol += np.cos((i+1)*x[i])
        return sol
        # return -sum(x)

    ncols = 3
    nrows = 4
    x_grid = np.linspace(0,1,ncols)
    y_grid = np.linspace(5,6,nrows)
    grid_x,grid_y = np.meshgrid(x_grid,y_grid)
    ntech = 20
    sol, arr = hybrid_sweep(obj_function,nrows,ncols,ntech)
    print(sol)
    print(arr)
    # x = np.linspace(0,1,10)
    # y = np.linspace(0,1,3)
    # xx,yy = np.meshgrid(x,y)
    # print(np.shape(xx))
    # print(xx)