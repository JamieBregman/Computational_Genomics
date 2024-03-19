#Week 4 Part 2 Assignment
#Jamie Bregman
#Python program to calculate N? (where "?" is a value betwen 10 - 90 that assesses the contiguity of a genome assembley)

#Importing modules
import math

filename = input("Please enter name of fasta file of assembled contigs: ")
N_value = int(input("Please enter N? value, where N is a value between 10 - 90: "))
genome_size = int(input("Please enter actual size of the genome in base pairs: "))

#Reading in fasta file of assembled contigs
def fasta_to_dict(filename):
    dictionary = {}
    current_species = ''
    current_sequence = []
    for line in open(filename):
        if line.startswith(">") and current_species == '':
            current_species = line.split(' ')[0]
        elif line.startswith(">") and current_species != '':
            dictionary[current_species] = ''.join(current_sequence)
            current_species = line.split(' ')[0]
            current_sequence = []
        else:
            current_sequence.append(line.rstrip())
    dictionary[current_species] = ''.join(current_sequence)
    return dictionary

#Calling function to open fasta file
my_dict = fasta_to_dict(filename)

#Creating a list from the sequences
contig_list = list(my_dict.values())

#Sorting contigs into a list of contigs from longest to shortest
contig_list.sort(key=len)
contig_list_sorted = contig_list[::-1]

#Creating one combined sequence list from all of the contigs
contig_list_sorted_and_joined = ["".join(contig_list_sorted)]

#Creating a dictionary where the keys are the indexes for the individual contigs in the contig list (contig_list)
#The values for each key are then the associated indexes within the combined list (contig_list_sorted_and_joined)
#For example, if the longest contig contained 1000 bases and the second largest contig contained 500 bases, then the dictionary would look like this:
#{0:[0,1,2...999], 1:[1000,1001,1002...1499], etc}
index_dict = {}
start = 0
for i, j in enumerate(contig_list_sorted):
    index = contig_list_sorted_and_joined[0].find(j, start)
    end_index = index + len(j)
    index_dict[i] = list(range(index, end_index))
    start = end_index

#Finding the index for the N? value based on the user inputs
#The exact base in the combined sequence we want is the base that pushes the length of contigs over the desired N value
#This base represents the index of the desired contig
def find_N_index(N_value,genome_size):
    
    #Find the number of bases until N value is reached (e.g., if N = 50, then N_index = the halfway point of the entire sequence)
    N_index = float((N_value/100) * genome_size)

    #If index is an integer, add one more to the index to find the base that pushes it over
    if (N_index).is_integer():
        N_index = N_index + 1
    #If index is not an integer, round it up for the same reason as above
    else:
        N_index = math.ceil(N_index)

    return N_index

#Calling the above function with the user inputted values
index = find_N_index(N_value,genome_size)

#Using the index found above and the dictionary created earlier with each contig's indexes in the full sequence to identify the contig that contains the desired base
key = [k for k,v in index_dict.items() if index in v]

contig_name = list(my_dict.keys())[list(my_dict.values()).index(contig_list[key])]

#Writing results to an output file
output_file = open("N?_output.txt","w")
#Writing the name of the contig and its length
output_file.write("The N"+str(N_value)+" value for "+str(filename)+" is:"+"\n" +
                  "\n" + list(my_dict.keys())[list(my_dict.values()).index(contig_list_sorted[int(str(key)[1:-1])])] + str(len(contig_list_sorted[int(str(key)[1:-1])])))
                  
output_file.close()
                  
