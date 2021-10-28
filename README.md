# Deep_Mutational_Scanning

This repository includes two basic functions:
- Deep Mutational Scanning: Generate a single-point mutation with defined amino acids upon a given sequence. Here the function AA_mutagenesis() takes:
    name_query;
    sequence_query of the fasta entry
    
- Heatmap display:  Generate a heatmap to display the results e.g. predictions applied to the library of single mutants from deep mutational scanning (sequence_query (seq_0) mutants):
    sequence_query: string of the reference peptide;
    df : Large-scale mutational data = dataframe of HAPMOD results of sequence_query (seq_0) mutants;
    fit_func: column of interest in df with HAPMOD results we want to display;
    Boolean True or False - if True, calculate differences with fitness value of sequence_query
    vmin, vmax: minimal and maximal values of fit_func (Boolean True: default -1,1 and False: default 0,1) 
