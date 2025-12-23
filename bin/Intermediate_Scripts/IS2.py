#!/usr/bin/env python
import sys
from Bio import SeqIO

#defining load_fasta function
def load_fasta(fasta_file):
    sequences = {} #start empty
    #reads each record in the fasta file using SeqIO module
    for record in SeqIO.parse(fasta_file, "fasta"):
    #save the sequences as the values, and the IDs as the keys
        sequences[record.id] = str(record.seq)
    return sequences

#defining comparison function
def compare_sequences(fasta1, fasta2, output_file):
    sequences1 = load_fasta(fasta1) #load transdecoderpep file
    sequences2 = load_fasta(fasta2) #load maturefasta

    with open(output_file, "w") as out_f:
        for seq_id, seq1 in sequences1.items():
            if seq_id in sequences2:
                seq2 = sequences2[seq_id]
                #see if the mature sequence from the mature fasta matches the end part of any full sequence in the transdecoder file
                if seq1.endswith(seq2):
                    diff = seq1[:-len(seq2)] #if there is a match then save the "difference" i.e. the signal sequence as diff
                else:
                    diff = ""  #no suffix match then leave it out

                if diff:  # Difference has to be not empty.
                    out_f.write(f">{seq_id}\n{diff}\n") #fasta format, > seq id new line, sequence

#makes sure we have 4 arguments
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python compare_fasta.py <input1.fasta> <input2.fasta> <output.fasta>")
        sys.exit(1)
        
#defines what arugment is what
    fasta_file1 = sys.argv[1]
    fasta_file2 = sys.argv[2]
    output_file = sys.argv[3]
#run the command 
    compare_sequences(fasta_file1, fasta_file2, output_file)

