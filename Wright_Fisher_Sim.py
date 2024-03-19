#Week 6 Assignment
#Wright-Fisher simulator with a diploid population that is hermaphroditic but obligately outcrossing. Includes discrete generations with a constant population
#Track for 2 different populations starting from the same point


#Importing modules
import random
import numpy

#Inputting user data (population size, number of generations to simulate, and an output file to write the simulation to)
population_size = int(input("Please enter population size: "))
num_generations = int(input("Please enter a number of generations to simulate: "))
output_file = input("Please enter name of output file: ")

#Function that creates an ongoing list of "mutation alleles"
#Essentially a running list of numbers that are introduced as new mutations into the population
#Exists so that once a mutation occurs, if it is lost to genetic drift than it does not reappear (ASSUMPTION)
def mutation_alleles_list():
    number = 1
    while True:
        yield number
        number +=1

#Calling the above function to generate the running list of mutation alleles
mutation_alleles = mutation_alleles_list()

#Function that creates a parental generation
#Takes the user inputted population size (n) as the input
#Generates a population of size n where each individual has a randomly generated pair of 1's and 0's as their genome
def create_parent_gen(population_size):
    parent_gen = []
    parent_alleles = [0,1]
    for i in range(population_size):
        parent_gen.append([random.choice(parent_alleles),random.choice(parent_alleles)])
    return parent_gen

#Calling the above function to generate a parent generation
parent_gen = create_parent_gen(population_size)

#Function to simulate the offspring generation from the parent generation
#Each offspring generation has the same size as the parent (ASSUMPTION)
#For each offspring, randomly chooses 2 parents WITH replacement. For each of those parents, that offspring randomly inherits one allele from each parent.
#Additionally includes a mutation rate for each offspring.
#When an offspring is "generated", a number is randomly generated in a normal distribution.
#If that number is greater than 2.5 SD's from the mean (~1% chance of mutation), then that particular offspring has a new mutation.
#The result is a new generation of offspring that have mostly the same genomes as the parents, but with a slight chance of mutation
def simulate_next_gen(parent_gen):
    new_gen = []
    for i in range(len(parent_gen)):
        genome1,genome2 = random.sample(parent_gen, k=2)
        new_gen.append([random.choice(genome1),random.choice(genome2)])

    for i in range(len(new_gen)):
        for j in range(len(new_gen[i])):
            probability = abs(numpy.random.normal())
            if probability >= 2.5:
                new_gen[i][j] = next(mutation_alleles)
            else:
                continue
    return new_gen

#Function to count unique items, which is used to report allele frequencies in the output file.
#This information is useful if one wanted to plot the frequency vs generation for a particular allele
def count_unique_items(list_of_lists):
    flat_list = [item for sublist in list_of_lists for item in sublist]
    
    frequencies = {}
    for item in flat_list:
        if item not in frequencies:
            frequencies[item] = 1
        else:
            frequencies[item] += 1
    
    total_items = len(flat_list)
    proportions = {k: v / total_items for k, v in frequencies.items()}
    
    return proportions

#Final for loop that creates new offspring from the original parent generation (where each subsequent offspring generation BECOMES the parent generation) for a user-input number of generations
def wf(parent_gen, num_generations):

    f = open(output_file,"a")
    
    for p in range(2):
        parent_gen1 = parent_gen
        allele_dict = count_unique_items(parent_gen1)
        f.write("population " + str(p) + "\n" + "\n")
        f.write(str(parent_gen1) + "\n" + "allele frequencies: " + "\n")
        for key, value in allele_dict.items():
            f.write(f'{key}: {value}\n')
        for i in range(num_generations):
            parent_gen1 = simulate_next_gen(parent_gen1)
            allele_dict = count_unique_items(parent_gen1)
            f.write(str(parent_gen1) + "\n" + "allele frequencies: " + "\n")
            for key, value in allele_dict.items():
                f.write(f'{key}: {value}\n')
    f.close()

wf(parent_gen, num_generations)
