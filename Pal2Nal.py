#Week 3 Part 3
#Jamie Bregman

#Inputting alignment & sequence file
#Files must be in fasta format and in the same order
alignment_file = input("Please enter alignment file in fasta format (items must be in the same order as sequence file): ")
sequence_file = input("Please enter DNA or mRNA sequence file in fasta format (items must be in the same order as alignment file): ")
output_file = input("Please enter name of output file: ")

#Function to import fasta file(s) as dictionaries
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

#Creating a dictionary for both the alignment file and the sequence file
alignment_dictionary = fasta_to_dict(alignment_file)
sequence_dictionary = fasta_to_dict(sequence_file)

#Creating a list of each alignment in the alignment file
alignment_list = list(alignment_dictionary.values())

#Creating a list of each sequence in the sequence file
sequence_list = list(sequence_dictionary.values())

#Creating a list of sequence names from the sequence file
sequence_names = list(sequence_dictionary.keys())

#Creating an empty list to hold the sequence list as codons
sequence_list_codons = []
#For loop that breaks each sequence into groups of 3 (codons) and appends them to the empty list created above
for i in sequence_list:
    sequence_list_codons.append([i[j:j+3] for j in range(0,len(i),3)])

#Dictionary of all possible codon combinations for each amino acid letter
amino_acid_dictionary ={
    "F": ["TTT","TTC"],
    "L": ["TTA","TTG","CTT","CTC","CTA","CTG"],
    "I": ["ATT","ATC","ATA"],
    "M": ["ATG",],
    "V": ["GTT","GTC","GTA","GTG"],
    "S": ["TCT","TCC""TCA","TCG","AGT","AGC"],
    "P": ["CCT","CCC","CCA","CCG"],
    "T": ["ACT","ACC","ACA","ACG"],
    "A": ["GCT","GCC","GCA","GCG"],
    "Y": ["TAT","TAC"],
    "H": ["CAT","CAC"],
    "Q": ["CAA","CAG"],
    "N": ["AAT","AAC"],
    "K": ["AAA","AAG"],
    "D": ["GAT","GAC"],
    "E": ["GAA","GAG"],
    "C": ["TGT","TGC"],
    "W": ["TGG",],
    "R": ["CGT","CGC","CGA","CGG","AGA","AGG"],
    "G": ["GGT","GGC","GGA","GGG"]}

#Pal2Nal example data had numbers in their sequences, so these are added in case input files have numbers
numbers = ["1","2","3","4","5","6","7","8","9","0"]

#Main for loop
#For each sequence in the MSA
for i in range(len(alignment_list)):
    #For each index in that sequence
    for j in range(len(alignment_list[i])):
        #Replacing all non amino acid characters with the appropriate replacement of length 3.
        if alignment_list[i][j] == "-":
            sequence_list_codons[i].insert(j,"---")
        elif alignment_list[i][j] == "_":
            sequence_list_codons[i].insert(j,"___")
        elif alignment_list[i][j] == "X":
            sequence_list_codons[i].insert(j,"XXX")
        elif alignment_list[i][j] == "*":
            sequence_list_codons[i].insert(j,"TAG")
        elif alignment_list[i][j] in numbers:
            sequence_list_codons[i].insert(j,"___")
        #If the sequence has an amino acid letter
        else:
            #Make sure that letter corresponds to the correct amino acid
            if sequence_list_codons[i][j] in amino_acid_dictionary[alignment_list[i][j]]:
                continue
            #If that letter does not correspond to the correct amino acid, create an error message that will be written to the output file
            #Including this because the example on PAL2NAL's website has a similar feature
            else:
                error_messages = []
                error_messages.append("WARNING: " + str(sequence_names[i]) + "at codon position " + str(j) + ": " + str(alignment_list[i][j]) + " does not correspond  to " + str(sequence_list_codons[i][j]))
                continue

#Joining the codons in each sequence to create the final codon alignment
codon_aligned_sequences = []
for k in sequence_list_codons:
    codon_aligned_sequences.append("".join(k))

#Creating output file for the program
with open (output_file,"w") as f:
    #Writing each sequence name and the corresponding sequence
    for l in range(len(sequence_names)):
        f.write(str(sequence_names[l]) + "\n" + str(codon_aligned_sequences[l]) + "\n")
    f.write("\n")
    #Writing any potential error messages
    for m in error_messages:
        f.write(f"{m}\n")
    f.close()
