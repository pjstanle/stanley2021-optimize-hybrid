import numpy as np
from math import log

class GeneticAlgorithm():
    """a simple genetic algorithm"""

    def __init__(self):

        # inputs
        self.bits = np.array([]) # array of ints same length as design_variables. 
        self.bounds = np.array([]) # array of tuples same length as design_variables
        self.variable_type = np.array([]) # array of strings same length as design_variables ('int' or 'float')
        self.objective_function = None # takes design_variables as an input and outputs the objective values (needs to account for any constraints already)
        self.max_generation = 100
        self.population_size = 0
        self.crossover_rate = 0.1
        self.mutation_rate = 0.01
        self.tol = 1E-6
        self.convergence_iters = 5
        
        # internal variables, you could output some of this info if you wanted
        self.design_variables = np.array([]) # the desgin variables as they are passed into self.objective function
        self.nbits = 0 # the total number of bits in each chromosome
        self.nvars = 0 # the total number of design variables
        self.parent_population = np.array([]) # 2D array containing all of the parent individuals
        self.offspring_population = np.array([]) # 2D array containing all of the offspring individuals
        self.parent_fitness = np.array([]) # array containing all of the parent fitnesses
        self.offspring_fitness = np.array([]) # array containing all of the offspring fitnesses
        self.discretized_variables = {} # a dict of arrays containing all of the discretized design variable

        # outputs
        self.solution_history = np.array([])
        self.optimized_function_value = 0.0
        self.optimized_design_variables = np.array([])


    def initialize_population(self):

        self.parent_population = np.random.randint(0,high=2,size=(self.population_size,self.nbits))
        self.offspring_population = np.zeros_like(self.parent_population)


    def chromosome_2_variables(self,chromosome):  
        """convert the binary chromosomes to design variable values"""      

        first_bit = 0

        float_ind = 0

        for i in range(self.nvars):
            binary_value = 0
            for j in range(self.bits[i]):
                binary_value += chromosome[first_bit+j]*2**j
            first_bit += self.bits[i]

            if self.variable_type[i] == "float":
                self.design_variables[i] = self.discretized_variables["float_var%s"%float_ind][binary_value]
                float_ind += 1

            elif self.variable_type[i] == "int":
                self.design_variables[i] = self.bounds[i][0] + binary_value

    
    def optimize_ga(self,print_progress=True):
        """run the genetic algorithm"""

        # print("start GA")
        # determine the number of design variables and initialize
        self.nvars = len(self.variable_type)
        self.design_variables = np.zeros(self.nvars)
        float_ind = 0
        for i in range(self.nvars):
            if self.variable_type[i] == "float":
                ndiscretizations = 2**self.bits[i]
                self.discretized_variables["float_var%s"%float_ind] = np.linspace(self.bounds[i][0],self.bounds[i][1],ndiscretizations)
                float_ind += 1

        # determine the total number of bits
        for i in range(self.nvars):
            if self.variable_type[i] == "int":
                int_range = self.bounds[i][1] - self.bounds[i][0]
                int_bits = int(np.ceil(log(int_range,2)))
                self.bits[i] = int_bits
            self.nbits += self.bits[i]        

        # initialize the population
        # print("initialize population")
        if self.population_size%2 == 1:
            self.population_size += 1

        self.initialize_population()

        # initialize the fitness arrays
        # print("initialize fitness")
        self.parent_fitness = np.zeros(self.population_size)
        self.offspring_fitness = np.zeros(self.population_size)

        # initialize fitness of the parent population
        for i in range(self.population_size):
            self.chromosome_2_variables(self.parent_population[i])
            self.parent_fitness[i] = self.objective_function(self.design_variables)

        converged = False
        ngens = 1
        generation = 1
        difference = self.tol * 10000.0
        self.solution_history = np.zeros(self.max_generation+1)
        self.solution_history[0] = np.min(self.parent_fitness)

        # print("start optimization")
        while converged==False and ngens < self.max_generation:
            self.crossover()
            self.mutate()
            for i in range(self.population_size):
                self.chromosome_2_variables(self.offspring_population[i])
                self.offspring_fitness[i] = self.objective_function(self.design_variables)

            # rank the total population from best to worst
            total_fitness = np.append(self.parent_fitness,self.offspring_fitness)
            ranked_fitness = np.argsort(total_fitness)[0:int(self.population_size)]

            # take the best. Might switch to some sort of tournament, need to read more about what is better
            # for now I've decided to only keep the best members of the population. I have a large population in 
            # the problems I've run with this so I assume sufficient diversity in the population is maintained from that
            total_population = np.vstack([self.parent_population,self.offspring_population])
            self.parent_population[:,:] = total_population[ranked_fitness,:]
            self.parent_fitness[:] = total_fitness[ranked_fitness]
            
            # store solution history and wrap up generation
            self.solution_history[generation] = np.min(self.parent_fitness)

            if generation > self.convergence_iters:
                difference = self.solution_history[generation-self.convergence_iters] - self.solution_history[generation]
            else:
                difference = 1000
            if abs(difference) <= self.tol:
                converged = True
            
            # shuffle up the order of the population
            shuffle_order = np.arange(1,self.population_size)
            np.random.shuffle(shuffle_order)
            shuffle_order = np.append([0],shuffle_order)
            self.parent_population = self.parent_population[shuffle_order]
            self.parent_fitness = self.parent_fitness[shuffle_order]
            if print_progress==True:
                print(self.parent_fitness[0])

            generation += 1
            ngens += 1

        # Assign final outputs
        self.solution_history = self.solution_history[0:ngens]
        self.optimized_function_value = np.min(self.parent_fitness)
        self.chromosome_2_variables(self.parent_population[np.argmin(self.parent_fitness)])
        self.optimized_design_variables = self.design_variables


    def crossover(self):
        # Random crossover

        # set offspring equal to parents
        self.offspring_population[:,:] = self.parent_population[:,:]

        # mate conscutive pairs of parents (0,1),(2,3), ...
        # The population is shuffled so this does not need to be randomized
        for i in range(int(self.population_size/2)):
            # trade bits in the offspring
            crossover_arr = np.random.rand(self.nbits)
            for j in range(self.nbits):
                if crossover_arr[j] < self.crossover_rate:
                    self.offspring_population[2*i][j], self.offspring_population[2*i+1][j] = self.offspring_population[2*i+1][j], self.offspring_population[2*i][j]

    
    def mutate(self):
        # Randomly mutate bits of each chromosome
        for i in range(int(self.population_size)):
            # mutate bits in the offspring
            mutate_arr = np.random.rand(self.nbits)
            for j in range(self.nbits):
                if mutate_arr[j] < self.mutation_rate:
                    self.offspring_population[i][j] = (self.offspring_population[i][j]+1)%2


if __name__=="__main__":

    def simple_obj(x):
        return x[0]+x[1]

    def rosenbrock_obj(x):
        return (1-x[0])**2 + 100.0*(x[1]-x[0]**2)**2

    def ackley_obj(x):
        p1 = -20.0*np.exp(-0.2*np.sqrt(0.5*(x[0]**2+x[1]**2)))
        p2 = np.exp(0.5*(np.cos(2.*np.pi*x[0]) + np.cos(2.0*np.pi*x[1]))) + np.e + 20.0
        return p1-p2

    def rastrigin_obj(x):
        A = 10.0
        n = len(x)
        tot = 0
        for i in range(n):
            tot += x[i]**2 - A*np.cos(2.0*np.pi*x[i])
        return A*n + tot


    import matplotlib.pyplot as plt

    ga = GeneticAlgorithm()
    ga.bits = np.array([8,8])
    ga.bounds = np.array([(0.0,1.),(0.,1.2)])
    # ga.variable_type = np.array(["int","int"])
    ga.variable_type = np.array(["float","float"])
    ga.population_size = 20
    ga.max_generation = 100
    ga.objective_function = rastrigin_obj
    ga.crossover_rate = 0.1
    ga.mutation_rate = 0.01
    ga.convergence_iters = 25
    ga.tol = 1E-8

    ga.optimize_ga()
    print("optimal function value: ", ga.optimized_function_value)
    print("optimal design variables: ", ga.optimized_design_variables)
    print("nbits: ", ga.nbits)
    plt.figure(1)
    plt.plot(ga.solution_history)


    from mpl_toolkits.mplot3d import Axes3D
    X = np.arange(-5, 5, 0.02)
    Y = np.arange(-5, 5, 0.02)
    X, Y = np.meshgrid(X, Y)
    Z = np.zeros_like(X)
    for i in range(np.shape(Z)[0]):
        for j in range(np.shape(Z)[1]):
            Z[i][j] = rastrigin_obj(np.array([X[i][j],Y[i][j]]))
    
    # Plot the surface.
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X, Y, Z,linewidth=0, antialiased=False)

    plt.show()
    
