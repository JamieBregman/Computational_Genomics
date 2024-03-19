#Week 2 Part 1
#Jamie Bregman

#Python script to generate a homologous chromosome with simulated mutations, insertions, deletions, translocations, & inversions
#Random normal distribution is used for randomness (see code comments below for individual reasoning)

#Importing modules
import numpy
import random

#Getting chromosome filename and output filename
filename = input("Please enter fasta file for your chromosome: ")
output_file = input("Please enter name of OUTPUT fasta file: ")

#Opening fasta file and parsing into a sequence name and a sequence
with open(filename) as f:
    f_split = f.read().splitlines()
    seq = []
    for line in f_split:
        if line.startswith(">"):
            sequence_name = line
        else:
            seq.append(line)

sequence = "".join([str(i) for i in seq])

#Editing sequence name
sequence_name = sequence_name + " SIMULATED"

#Creating empty "map" for positions (indexes) of mutations, insertions, and deletions
mutation_or_indel_map = []

#Creating empty "map" for positions (indexes) of translocations and inversions
translocation_or_inversion_map = []

#Function to build the mutation maps
def get_mutation_index(sequence):

    #Looping over each base in the sequence
    for i in range(len(sequence)):

        #Generate a random number from a normal distribution (mean = 0, SD = 1)
        probability = abs(numpy.random.normal())

        #If the random number is between 2 - 3 SD's away from the mean, mark that index as a potential mutation/indel (> 5% chance)
        #Reasoning: low chance for mutation/indel, but still more common than a translocation or inversion
        if probability >= 2 and probability <= 3:
            mutation_or_indel_map.append(i)

        #If the random number is farther than 3 SD's away from the mean, mark that index as a potential translocation/inversion (> 0.5% chance)
        #Reasoning: very low chance for a translocation or inversion
        elif probability > 3:
            translocation_or_inversion_map.append(i)

        #All other indexes are marked as no mutation 
        else:
            continue

#Getting the mutation maps
get_mutation_index(sequence)

#List of bases that are used to create mutations later
bases = ["A","C","G","T"]

#Function that generates a random sequence of bases
#Used for insertions to insert random base(s)
#Used for translocations to simulate the addition of genetic material from another chromosome
def generate_random_sequence(length):
    return "".join(random.choices(bases, k = length))

#Function that creates a new mutated sequence based on the mutation maps and the original sequence
def create_new_sequence(sequence, mutation_or_indel_map, translocation_or_inversion_map):

    #Turning the original sequence into a "list" of bases so that they can be changed
    new_sequence = list(sequence)

    #For each index in the sequence...
    #Working backwards through the indexes to prevent index errors from mutations
    for i in reversed(range(len(sequence))):

        if i not in mutation_or_indel_map and i not in translocation_or_inversion_map:
            continue
        
        #If that index matches one of the indexes for a mutation or indel:
        elif i in mutation_or_indel_map:

            #Randomly choose type of mutation
            choices1 = ["mutation", "insertion", "deletion"]

            mutation_or_indel_choice = numpy.random.choice(choices1)

            #If mutation...
            if mutation_or_indel_choice == "mutation":
                #Randomly replace that base with a different base
                new_sequence[i] = numpy.random.choice([j for j in bases if j != sequence[i]])
                  
            #If insertion...
            elif mutation_or_indel_choice == "insertion":
                #Use random normal distribution (mean = 1, SD = 10) to simulate length of insertion
                #Reasoning: most insertions are small but the distribution offers some variation for potential larger insertions
                insertion_length = round(abs(numpy.random.normal(1,10)))
                    
                #Inserting the base(s) into the new sequence
                new_sequence.insert(i,generate_random_sequence(insertion_length))

                #If deletion...        
            elif mutation_or_indel_choice == "deletion":
                #Use same random normal distribution paramaters as above to simulate the length of the deletion (same reasoning as above)
                def calculate_deletion_length(sequence):
                    deletion_length = round(abs(numpy.random.normal(1,10)))

                    #While loop to ensure that for indexes at the end of the sequence, the length of the mutation does not exceed the length of the sequence
                    while True:
                        if i+deletion_length < len(sequence):
                            return deletion_length
                        else:
                            deletion_length = round(abs(numpy.random.normal(1,10)))

                    return deletion_length

                deletion_length = calculate_deletion_length(sequence)
                            
                #Deleting the base(s) in the new sequence
                del new_sequence[i:i+deletion_length]
                    
        #If that index matches one of the indexes for a translocation or inversion:
        elif i in translocation_or_inversion_map:

            #Similar to above, randomly choose which type of mutation
            choices2 = ["translocation", "inversion"]

            translocation_or_inversion_choice = numpy.random.choice(choices2)

            #If translocation...
            if translocation_or_inversion_choice == "translocation":
                #Use random normal distribution (mean = 100, SD = 100) to simulate length of translocation
                #Reasoning: since translocations are entire pieces of chromosome that are swapped, I widened the distribution to create much longer sequences
                def calculate_translocation_length(sequence):
                    translocation_length = round(abs(numpy.random.normal(100,100)))

                    #While loop to ensure that for indexes at the end of the sequence, the length of the mutation does not exceed the length of the sequence
                    while True:
                        if i+translocation_length < len(sequence):
                            return translocation_length
                        else:
                            translocation_length = round(abs(numpy.random.normal(1,10)))

                    return translocation_length

                translocation_length = calculate_deletion_length(sequence)

                #Deleting the appropriate length of the sequence
                del new_sequence[i:i+translocation_length]
                #Inserting a random sequence of the appropriate length to simulate genetic material from another chromosome
                new_sequence.insert(i,generate_random_sequence(translocation_length))

            #If inversion...
            elif translocation_or_inversion_choice == "inversion":
                #Use same random normal distribution paramaters as above to simulate the length of the inversion (same reasoning as above)
                def calculate_inversion_length(sequence):
                    inversion_length = round(abs(numpy.random.normal(100,100)))

                    #While loop to ensure that for indexes at the end of the sequence, the length of the mutation does not exceed the length of the sequence
                    while True:
                        if i+inversion_length < len(sequence):
                            return inversion_length
                        else:
                            inversion_length = round(abs(numpy.random.normal(1,10)))

                    return inversion_length

                inversion_length = calculate_deletion_length(sequence)

                #New sequence is equal to everything before the index + the specific index that you are inverting + everything that comes after that inverted region
                new_sequence = new_sequence[:i] + new_sequence[i:i+inversion_length][::-1] + new_sequence[i+inversion_length:len(new_sequence)]

    #Joining the various mutation and sequence strings into one continuous, simulated sequence
    return "".join(new_sequence)

#Calling the above function to create the final sequence
final_sequence = create_new_sequence(sequence, mutation_or_indel_map, translocation_or_inversion_map)

#Writing the final sequence to a file in fasta format with lines that are 70 base pairs long
with open(output_file, "w") as w:
    w.write(sequence_name + "\n")
    for i in range(0,len(final_sequence), 70):
        w.write(final_sequence[i:i+70] + "\n")

#Closing the output file
w.close()
