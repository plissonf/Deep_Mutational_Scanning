import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re

np.random.seed(19680801)

def AA_mutagenesis(name_query, sequence_query):
    '''Generate a single-point mutation with defined amino acids upon a given sequence. Here the function AA_mutagenesis()
    should take name_query and sequence_query of the fasta entry'''

    # Create two lists for new mutants ids and sequences, include sequence candidate (name_0, seq_0)
    name_0 = name_query
    seq_0 = sequence_query
    ids =[name_0,]
    seqs =[seq_0,]
    list_aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']


    def AA_mutation(sequence, position, mutation):
    # Generate a single-point mutation upon a given sequence
        mutation = str(mutation)
        for char in sequence:
            if char == position.upper() or char == position.lower():
                sequence = sequence.replace(char, mutation)
                return sequence

    # AA mutagenesis applied to seq_0, new ids with name_0
    for aa in list_aa:
        for i in range(len(seq_0)):
            mutant = AA_mutation(seq_0, seq_0[i], aa)
            new_name = '[' + seq_0[i] + str(i+1) + aa + ']' + name_0
            ids.append(new_name)
            seqs.append(mutant)
            #return ids, seqs
            #print(new_name, mutant)

    # Convert that labelled list into a dataframe
    mutagenesis_df = pd.DataFrame(seqs, columns=['SEQUENCE'], index=ids)
    link = "./Data/" + name_0 + '_mutagenesis.csv'
    mutagenesis_df.to_csv(link)

    def fasta_converter(database):
    # Convert that labelled list into a dataframe
        '''Use the indices of the database and the column 'SEQUENCE' to create a fasta file'''
        link2 = "./Data/" + name_0 + '_mutagenesis.fasta'
        ofile = open(link2, "w")

        for i in range(len(database.SEQUENCE)):
            ofile.write(">" + database.index[i] + "\n" +database.SEQUENCE[i] + "\n")

        ofile.close()

    # Convert that dataframe into fasta file
    fasta_converter(mutagenesis_df)

if __name__ == "__main__":
    AA_mutagenesis() 
