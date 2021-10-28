import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
import plotly as py
import simplejson as json

np.random.seed(19680801)

def Heatmap_mutagenesis(sequence_query, df, fit_func, diff, vmin, vmax):
    ''' Generate a heatmap to display the results e.g. predictions applied to the library of single mutants from deep mutational      scanning (sequence_query (seq_0) mutants):
    # sequence_query: string of the reference peptide;
    # df : Large-scale mutational data = dataframe of HAPMOD results of sequence_query (seq_0) mutants;
    # fit_func: column of interest in df with HAPMOD results we want to display;
    # Boolean True or False - if True, calculate differences with fitness value of sequence_query
    # vmin, vmax: minimal and maximal values of fit_func (Boolean True: default -1,1 and False: default 0,1) '''

    seq_0 = sequence_query
    #list of amino acids
    list_aa = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

    # Hide the first row (reference peptide)
    sub_df = df.iloc[1:]


    # Create a column 'AA point-mutation (AAPM)' that displays the nature of amino acid susbtituting the AA in a given position
    x_values = list_aa * len(seq_0)
    x_values.sort()
    sub_df['x'] = x_values

    # Create a column 'Original AA (OAA)' that displays the nature of amino acid susbtituting the AA in a given position
    list_oa = list(seq_0)
    y_values = list_oa * len(list_aa)
    sub_df['y'] = y_values

    # Subset dataframe based on model values of choice
    df_xyz = sub_df[['x', 'y', fit_func]]

    df_xyz['z_diffs'] = df_xyz[fit_func]

    if diff == True:
        ref_value = df[fit_func].iloc[[0]]
        df_xyz['z_diffs'] = df_xyz[fit_func] - float(ref_value)

    else:
        pass

    # Use pivot method to re-arrange the dataset, reorder columns based on reference sequence
    table = pd.pivot_table(df_xyz, index='x', columns='y', values='z_diffs')
    table = table[list_oa]

    # Measure mean values
    mean_rows = pd.DataFrame(table.mean(axis=1))
    mean_cols = pd.DataFrame(table.mean(axis=0))
    mean_cols = mean_cols.transpose()

    ## Plot heatmap
    fig = plt.gcf()
    fig.set_size_inches(10,8)
    grid = plt.GridSpec(nrows=2, ncols=3, width_ratios=[len(seq_0),1, 1], height_ratios=[len(list_aa), 1])


    # Main heatmap (ax1)
    ax1 = plt.subplot(grid[0])
    ax1 = sns.heatmap(table, linewidths=1, linecolor='Grey', cbar=False, cmap='Spectral_r', vmin=vmin, vmax=vmax)
    ax1.yaxis.tick_left()
    ax1.xaxis.tick_top()

    plt.ylabel("Amino Acid subtitution ", size=12)
    plt.yticks(rotation=0)
    plt.xlabel(" ", size=12)


    # Mean rows values
    ax2 = plt.subplot(grid[1])
    ax2 = sns.heatmap(mean_rows, linewidths=1, linecolor='Grey', xticklabels=False, cbar=False, cmap='Spectral_r', vmin=vmin, vmax=vmax)
    ax2.yaxis.tick_left()
    ax2.set_xlabel('$<ΔE(X)_i^x>_i$', labelpad=-400, fontsize=12, fontfamily='monospace')

    plt.ylabel("  ", size=12)
    plt.yticks(rotation=0)


    # Mean columns values
    ax4 = plt.subplot(grid[3])
    ax4 = sns.heatmap(mean_cols, linewidths=1, linecolor='Grey', yticklabels=False, cbar=False, cmap='Spectral_r', vmin=vmin, vmax=vmax)
    ax4.xaxis.tick_top()
    ax4.set_ylabel('$<ΔE(X)_i^x>_x$', labelpad=-430, rotation=0, fontsize=12, loc='center', fontfamily='monospace')

    plt.xlabel("Original Sequence", size=12)

    #Legend
    im = ax1.imshow(np.random.random((10,10)), cmap='Spectral_r', vmin=vmin, vmax=vmax)
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.78, 0.3, 0.02, 0.5])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(label='$ΔE(X)_i^x$', fontsize=12, rotation=270, fontfamily='monospace', loc='center', labelpad=20)

    #return plt.show()
    graphJSON = json.dumps(fig, cls=py.utils.PlotlyJSONEncoder)
    graphs.append(graphJSON)
    return graphs

if __name__ == "__main__":
    Heatmap_mutagenesis()

## Examples

#Effects of single-point mutagenesis of B2RP on HemoPI-1 model probabilities (differences)
#Heatmap_mutagenesis('GIWDTIKSMGKVFAGKILQNL', B2RP_mutants_models, 'md1.avg', True, -0.6, 0.6)

#Effects of single-point mutagenesis of B2RP on HemoPI-2 model probabilities (raw data)
#Heatmap_mutagenesis('GIWDTIKSMGKVFAGKILQNL', B2RP_mutants_models, 'md2.avg', False, 0, 1)

#Effects of single-point mutagenesis of B2RP on HemoPI-1 outlier scores (differences)
#Heatmap_mutagenesis('GIWDTIKSMGKVFAGKILQNL', B2RP_mutants_models, 'md1.OS', True, -0.4, 0.4)
