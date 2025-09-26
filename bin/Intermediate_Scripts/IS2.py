#!/usr/bin/env python
import sys
from Bio import SeqIO

#defining our fasta file loading function
def load_fasta(fasta_file):
    sequences = {} #start empty
    for record in SeqIO.parse(fasta_file, "fasta"): #reads each record in the fasta file using SeqIO module
    #save the sequences as the values, and the IDs as the keys
        sequences[record.id] = str(record.seq)
    return sequences

#defining our function to compare the sequences with the same header between two files and save their "difference" as an output file
def compare_sequences(fasta1, fasta2, output_file):
    sequences1 = load_fasta(fasta1) #loading in our transdecoder pep file
    sequences2 = load_fasta(fasta2) #loading in our mature.fasta file

    with open(output_file, "w") as out_f:
        for seq_id, seq1 in sequences1.items():
            if seq_id in sequences2:
                seq2 = sequences2[seq_id]

                # Find prefix or suffix difference
                if seq1.endswith(seq2):  #this checks if the mature sequence that matches the end of the sequence in the full file
                    diff = seq1[:-len(seq2)] #now we save everything in that transdecoder peptide file right before what's left would be the mature sequence the length of the mature sequence.
                else:
                    diff = ""  #no suffix match then leave it out

                if diff:  # Difference has to be not empty.
                    out_f.write(f">{seq_id}\n{diff}\n") #fasta format, > seq id new line, sequence

#makes sure we have 4 arguments
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python compare_fasta.py <input1.fasta> <input2.fasta> <output.fasta>")
        sys.exit(1)
        
#defines what arugment is whattttt
    fasta_file1 = sys.argv[1]
    fasta_file2 = sys.argv[2]
    output_file = sys.argv[3]
#run the command 
    compare_sequences(fasta_file1, fasta_file2, output_file)

