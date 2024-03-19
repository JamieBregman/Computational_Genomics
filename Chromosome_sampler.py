#Week 3 Part 2
#Jamie Bregman

#Importing modules
import random
import math

#Obtaining program arguments from command line
input_file = input("Please enter name of fasta input file: ")
output_file = input("Please enter name of fastq output file: ")
length_of_reads = int(input("Please enter length of sample reads: "))
coverage = int(input("Please enter coverage level (e.g. 10): "))

#Main function (with 4 arguments)
def main(input_file, output_file, length_of_reads, coverage):

    #Empty list to store fasta sequence
    sequence = []

    #Opening the fasta file
    with open(input_file) as f:
        f_split = f.read().splitlines()
        for line in f_split:
            #Storing the sequence name a variable
            #Storing everything else (the sequence) as a separate variable
            if line.startswith(">"):
                sequence_name = line[:11]
            else:
                sequence.append(line)

    #Editing sequence name for fastq format
    sequence_name_final = sequence_name.replace(">","@")

    #Joining the lines of the fasta file into one long sequence (strand1)
    strand1 = "".join([str(i) for i in sequence])

    #Creating a reverse complement strand (strand2)
    reverse = strand1[::-1]
    NoG=reverse.replace("G","X")
    NoC=NoG.replace("C","G")
    NoX=NoC.replace("X","C")
    NoT=NoX.replace("T","Q")
    NoA=NoT.replace("A","T")
    strand2 = NoA.replace("Q","A")

    #Function that pulls out a random substring from a larger substring
    #Length of substring is equal to the length of reads argument
    #Allows you to randomly sample the entire sequence for smaller reads
    def random_substring(seq):
        length = len(seq)
        start = random.randint(0,length - length_of_reads)
        return seq[start:start+length_of_reads]

    #Creating an empty list to store these samples of length 100 bp
    samples = []

    #Function to actually create the strands
    #According to the equation for read depth (R = YZ/X), the amount of reads necessary can be found with the equation Z = RX/Y.
    #The number of reads (Z) is equal to the (coverage (R) * length of the sequence (X)) / (length of samples (Y))
    #In this case, the coverage is 10X, the length of the sequence is 100,000 bp, and the length of reads is Y, giving us (10*10000)/(100) = 10,000
    #Because I wanted to sample evenly from both strands, I created a while loop that randomly samples (using the random_substring method above) 5,000 times from each strand for a total of 10,000 total reads.
    def create_samples(strand):
        n = 0
        while n < ((coverage*len(strand))/length_of_reads)/2:
            samples.append(random_substring(strand))
            n +=1

        return samples
    
    #Calling the method to create a long list of samples (10,000 samples) from both strands
    create_samples(strand1)
    create_samples(strand2)

    #Function that mutates the DNA to introduce an "error rate" in this sampling process
    #For each base in each sample (all 10,000), this function randomly generates a number between 0 - 1. If that probability is below 0.0055, then that base is randomly changed.
    #The reason 0.0055 is used is because this literature suggested this was the average error rate per base for Illumina sequencing (https://www.pnas.org/doi/10.1073/pnas.1319590110)
    def mutate_dna(dna):
        base_pairs = "ACGT"

        for i in dna:
            for j in i:
                if random.randint(0,1) <= 0.0055:
                    i.replace(j,random.choice(base_pairs))
                else:
                    continue

        return dna

    #Calling the error rate function
    #Currently have a list of 10000 reads that are 100 bp long and have errors ata rate of 0.0055 per base pair
    mutate_dna(samples)

    #Function to assign quality scores
    #Creates an empty quality score string that gets created as you iterate over a read
    #Used a random uniform distribution to assign Phred scores, in which the random number generates an error probability that gets placed into the Phred score equation to generate the appropriate character
    def assign_quality_scores(dna):
        quality_scores = ""
        for i in dna:
            error_prob = random.uniform(0,1)
            phred_score = -10 * (math.log10(error_prob))
            score = chr(int(round(phred_score)) + 33)
            quality_scores += score
        return quality_scores

    #Creating the output file
    output = open(output_file,"w")

    #For loop that loops over each read that has been created
    #Write the following (fastq format):
    #sequence name
    #sequence
    #+
    #quality scores
        
    for i in samples:
        output.write(sequence_name_final+"\n"+
                     str(i)+"\n"+
                     "+"+"\n"+
                     str(assign_quality_scores(i))+"\n")

    #Closing output file
    output.close()

#Calling main function
main(input_file, output_file, length_of_reads, coverage)
