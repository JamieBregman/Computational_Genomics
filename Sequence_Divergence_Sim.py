#Week 7 Assignment
#Simulating sequence divergence in a population with mutational effects on gene expression

#Importing modules
import random
import numpy
import math
import re

#User input (population size (constant for every generation), number of generations to simulate, and output file name)
population_size = int(input("Please enter population size: "))
num_generations = int(input("Please enter a number of generations to simulate: "))
output_file = input("Please enter name of output file: ")

#Function that creates an "individual" by creating a basic genome structure with promoter, enhancer, and repressor regions
def create_genome():
    num_enhancer = math.ceil(abs(numpy.random.normal())) #Random normal distribution models the number of enhancers (rounded up to avoid 0's)
    enhancer_vals = [] #Empty list to store enhancer values
    for i in range(num_enhancer): #For each enhancer
        enhancer_val = round(random.uniform(1,9)) #Randomly assign a score from 1 - 9
        #The below 4 lines of code add the enhancers to a final enhancer value list in a standardized fashion (ex: E3 for enhancer with effect 3)
        enhancer_vals.append("E")
        enhancer_vals.append(enhancer_val)
    enhancer_vals = [str(j) for j in enhancer_vals]
    enhancer_vals = "".join(enhancer_vals)

    #Promoter is modeled exactly the same as enhancer (above)
    num_promoter = math.ceil(abs(numpy.random.normal()))
    promoter_vals = []
    for i in range(num_promoter):
        promoter_val = round(random.uniform(1,9))
        promoter_vals.append("P")
        promoter_vals.append(promoter_val)
    promoter_vals = [str(j) for j in promoter_vals]
    promoter_vals = "".join(promoter_vals)

    #Repressor is modeled exactly like enhancer & promoter, except for 2 changes noted below
    num_repressor = round(abs(numpy.random.normal())) #Random normal distribution models the number of repressors (rounded up OR down to include 0's)
    repressor_vals = [] 
    for i in range(num_repressor):
        repressor_val = round(random.uniform(0.1, 0.99),2) #Repressor score is distributed as a decimal. This means when scores are multiplied, repressor will NEGATIVELY affect fitness
        repressor_vals.append("R")
        repressor_vals.append(repressor_val)
    repressor_vals = [str(j) for j in repressor_vals]
    repressor_vals = "".join(repressor_vals)

    #Gene is modeled exactly as the above 3 regulatory regions, except there is no model for how many genes there are. This stays constant at 1.
    gene_vals = []
    for i in range(1):
        gene_val = round(random.uniform(1,9))
        gene_vals.append("G")
        gene_vals.append(gene_val)
    gene_vals = [str(j) for j in gene_vals]
    gene_vals = "".join(gene_vals)

    #Genome is returned as a neat string that combines the above regulatory regions
    gene = enhancer_vals + promoter_vals + repressor_vals + gene_vals
    gene = "".join(gene)
    return gene
    
#Function that creates a parent generation of randomly chosen diploid genomes
def create_parent_gen(population_size):
    parent_gen = [] #Empty string to hold parent generation
    for i in range(population_size):
        individual = [create_genome(),create_genome()] #Create x new individuals with 2 randomly created genes, where x is equal to user-input population size
        parent_gen.append(individual)
        
    return parent_gen

parent_gen = create_parent_gen(population_size) #Calling the above function to create parent generation (generation 0)

#Function that splits apart a genome into its individual regions and returns them as a list of lists
def split_genome(genome):
    result = []
    i = 0
    while i < len(genome):
        letter = genome[i]
        i += 1
        
        number = ""
        while i < len(genome) and (genome[i].isdigit() or genome[i] == "."):
            number += genome[i]
            i += 1
        result.append([letter, number])
    return result

#Function that mutates a genome
def mutate_genome(individual):
    
    split = split_genome(individual) #Splitting apart the genome 
    for item in split: 
        if item[0] == "E" or item[0] == "P" or item[0] == "G": #For enhancers, promoters, and the gene
            mutation_value = abs(numpy.random.normal()) #Whether mutation is positive or negative

            if mutation_value >= 0 and mutation_value <= 1: #if mutation_value is between 0 to 1, mutation DECREASES (more common)
                mutation_change = round(random.uniform(1, float(item[1]))) #Randomly generate a value with which the region will change
                item[1] = int(item[1]) - mutation_change #initiate that change
                                
            elif mutation_value > 1: #if mutation_value is greater than 1, mutation INCREASES
                mutation_change = round(random.uniform(1,5))
                item[1] = int(item[1]) + mutation_change
                
        elif item[0] == "R": #For repressors
            mutation_value = abs(numpy.random.normal()) #Whether mutation is positive or negative

            if mutation_value >= 0 and mutation_value <= 1: #if mutation_value is between 0 to 1, mutation DECREASES (more common)
                mutation_change = round(random.uniform(0.1,float(item[1]))) 
                item[1] = float(item[1]) - mutation_change
                                
            elif mutation_value > 1: #if mutation_value is greater than 1, mutation INCREASES
                mutation_change = round(random.uniform(0.1, 0.99),2)
                item[1] = float(item[1]) + mutation_change

    #Removing any negative or elements that are 0, essentially an extra negative mutation
    for item in split:
        if float(item[1]) <= 0:
            del(item)
        else:
            continue

    #Formatting the final string
    mutated_genome = sum(split,[])
    mutated_genome_final = "".join(str(i) for i in mutated_genome)

    return mutated_genome_final

#Function that calculates organismal fitness from total genome expression
def get_fitness(individual):
    fitness_vals = re.findall(r'\d+\.*\d*', individual) #Uses re module to separate the letters from the numbers in a gene string to isolate only the numbers
    fitness_vals = [float(j) if '.' in j else int(j) for j in fitness_vals] #Turning the numbers back into floats and ints 

    score = numpy.prod(fitness_vals) #Multiplying every value in a gene string to calculate expression level for one individual gene
    return score

#Function that simulates the next generation
def simulate_next_gen(parent_gen):
    #Empty lists for new generation and the weights used to select parent organisms
    new_gen = [] 
    weights = []
    weights_final = []

    for i in parent_gen: #For every individual
        for j in i: #For both genes in that individual's genome
            fitness_score = get_fitness(j) #Get the expression level of that gene
            weights.append(fitness_score) #Append that level to the running list of "weights"
    weights = [weights[i] + weights[i+1] for i in range(0, len(weights)-1,2)] #Adding expression levels for both genes in an organism's genome to create one total expression level

    for i in weights: 
        weights_final.append(i/sum(weights)) #Creating a "final weights" list, where the expression level relative to every other organism is used as a weighting factor
                                             #In other words, the higher an expression level, the higher the fitness, and the more likely it is to be chosen as a parent for the next generation
    
    for i in range(len(parent_gen)):
        indexes = [i for i in range(len(parent_gen))] #Creating a list of indexes for the parent generation for easier choosing of parental genomes
        genome1,genome2 = numpy.random.choice(indexes, size=2, replace=False, p = weights_final) #Selecting WITH weight and WITHOUT replacement 2 parental genomes 
        new_gen.append([random.choice(parent_gen[genome1]),random.choice(parent_gen[genome2])]) #Randomly choose one gene from each parent chosen to pass down to the offspring

    for i in range(len(new_gen)):
        for j in range(len(new_gen[i])): #For every "offspring", model a mutation
            probability = abs(numpy.random.normal()) 
            if probability >= 2.5: #If probabilty from a random normal distribution is > 2.5, then mutate that offspring
                mutation = mutate_genome(new_gen[i][j])
                new_gen[i][j] = mutation #Mutation occurs by replacing one of the offspring's genes with a new gene created from paren
            else:
                continue

    return new_gen

#Function that calculates distribution of expression levels in a given generation
#Essentially the same code used to calculate the weights as seen above
#Expression is calculated by multiplying each regulatory region in a gene, and then adding both genes from a genome together to get total organismal expression
def calc_expression_levels(generation):
    expression_levels = []
    for i in generation:
        for j in i:
            fitness_score = get_fitness(j)
            expression_levels.append(fitness_score)
    expression_levels = [expression_levels[i] + expression_levels[i+1] for i in range(0, len(expression_levels)-1,2)]

    return expression_levels

#Function that runs the simulator for a user-input number of generations and writes the simulation to an output file
def simulation(parent_gen, num_generations):
    f = open(output_file,"a") #Creating output file

    for p in range(2): #Runs on a range of 2 to simulate population divergence (can be changed to see divergence in any number of populations)
        parent_gen1 = parent_gen #Making starting generation the same for both populations
        expression_distribution = calc_expression_levels(parent_gen1)
        f.write("population " + str(p) + "\n" + "\n")
        f.write(str(parent_gen1) + "\n" + "distribution of expression levels: " + str(expression_distribution) + "\n")
        for x in range(num_generations): #For user-input number of generations
            parent_gen1 = simulate_next_gen(parent_gen1) #keep simulating new generations
            expression_distribution = calc_expression_levels(parent_gen1) #Calculate expression distribution for each generation
            f.write(str(parent_gen1) + "\n" + "distribution of expression levels: " + str(expression_distribution) + "\n") #Write both the generation and its expression distribution to the file
    f.close()

simulation(parent_gen, num_generations) #Calling the main function
    
