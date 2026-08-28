'''
Functions to create the individual sequence output files. This
includes also the tabular .txt output of all positions.

This is implemented using multiprocessing, as matplotlib plotting
is the bottleneck for larger jobs.
'''
import numpy as np
import matplotlib
matplotlib.use('Agg') #this prevents failure during multiprocessing when matplotlib chooses another default.
matplotlib.set_loglevel('error') # this prevents PostScript backend warnings when using `eps`.
import matplotlib.pyplot as plt
from typing import List
from tqdm import tqdm
import multiprocessing
import os
import re


#convert label ids to simplified labels for plotting - only want one letter each
IDS_TO_PLOT_LABEL = {
                        0:  'O',
                        3:  'N', #SP
                        4:  'H',
                        5:  'C',
                        9:  'N', #LIPO
                        10: 'H',
                        12: 'c',
                        16: 'N', #TAT
                        17: 'R',
                        18: 'H',
                        19: 'C',
                        23: 'N', #TATLIPO
                        24: 'R',
                        25: 'H',
                        27: 'c',
                        31: 'P', #PILIN
                        32: 'c',
                        33: 'H',

}
#colors for letters, depending on region
REGION_PLOT_COLORS= {'N': 'red',
                     'H': 'orange',
                     'C': 'gold',
                     'I': 'gray',
                     'M': 'gray',
                     'O': 'gray',
                     'c': 'cyan',
                     'R': 'lime',
                     'P': 'red'
                     }
#colors for the probability lines, should somewhat match the letters
PROB_PLOT_COLORS = {3:  'red',
                    4:  'orange',
                    5:  'gold',
                    9:  'red',
                    10: 'orange',
                    12: 'cyan',
                    16: 'red',
                    17: 'lime',
                    18: 'orange',
                    19: 'gold',
                    23: 'red',
                    24: 'lime',
                    25: 'orange',
                    27: 'cyan',
                    31: 'red',
                    32: 'gold',
                    33: 'orange',

                }
# labels to write to legend for probability channels
IDX_LABEL_MAP = {
                 3:  'Sec/SPI n',
                 4:  'Sec/SPI h',
                 5:  'Sec/SPI c',
                 6:  'Sec/SPI I',
                 7:  'Sec/SPI M',
                 8:  'Sec/SPI O',
                 9:  'Sec/SPII n',
                 10: 'Sec/SPII h',
                 12: 'Sec/SPII cys',
                 13: 'Sec/SPII I',
                 14: 'Sec/SPII M',
                 15: 'Sec/SPII O',
                 16: 'Tat/SPI n',
                 17: 'Tat/SPI RR',
                 18: 'Tat/SPI h',
                 19: 'Tat/SPI c',
                 20: 'Tat/SPI I',
                 21: 'Tat/SPI M',
                 22: 'Tat/SPI O',
                 23: 'Tat/SPII n',
                 24: 'Tat/SPII RR',
                 25: 'Tat/SPII h',
                 27: 'Tat/SPII cys',
                 28: 'Tat/SPII I',
                 29: 'Tat/SPII M',
                 30: 'Tat/SPII O',
                 31: 'Sec/SPIII P', #PILIN
                 32: 'Sec/SPIII cons.',
                 33: 'Sec/SPIII h',
                 34: 'Sec/SPIII I',
                 35: 'Sec/SPIII M',
                 36: 'Sec/SPIII O',
}

#sum label ids into one other label, iterate over rest for plotting
LABEL_IDS_TO_PLOT = [3,4,5,9,10,12,16,17,18,19,23,24,25,27,31,32,33]


def sequence_plot(marginal_probs, viterbi_path, amino_acid_sequence, cleavage_site=-1, hide_threshold = 0.01, figsize = (12,4.5), title=None):
    
    amino_acid_sequence = amino_acid_sequence[:70]
    seq_length = len(amino_acid_sequence)
    pos_labels_to_plot = [IDS_TO_PLOT_LABEL[x] for x in viterbi_path[:seq_length]]
    joined_ticks = [x+'\n'+y for x,y in zip(pos_labels_to_plot, amino_acid_sequence)]
    tick_colors = [REGION_PLOT_COLORS[x] for x in pos_labels_to_plot]

    fig =plt.figure(figsize=figsize)

    for x in np.arange(1, seq_length+1):
        plt.axvline(x, c='whitesmoke') # Gridline at each pos
        tag = pos_labels_to_plot[x-1] # Predicted label at each pos
        plt.text(x,-0.05, tag, c=REGION_PLOT_COLORS[tag], ha='center', va='center')
        aa =  amino_acid_sequence[x-1] # AA at each pos
        plt.text(x,-0.1,aa, ha='center', va='center')


    # post-processed marginals: all other states summed in 0
    other_probs = marginal_probs[:seq_length,0]
    plt.plot(np.arange(1,seq_length+1), other_probs, label ='OTHER', c='lightcoral', linestyle='--')

    colors_used = [] #keep track of colors used - when already used, switch linestyle
    #NOTE this is still not optimal, only works when colors not used more than two times
    for idx in LABEL_IDS_TO_PLOT:
        probs = marginal_probs[:seq_length,idx]
        #skip low
        if (probs>hide_threshold).any():
            linestyle = '--' if PROB_PLOT_COLORS[idx] in colors_used else '-'
            plt.plot(np.arange(1,seq_length+1), probs, label = IDX_LABEL_MAP[idx], c=PROB_PLOT_COLORS[idx], linestyle=linestyle)
            colors_used.append(PROB_PLOT_COLORS[idx])


    if cleavage_site >0:
        plt.plot((cleavage_site, cleavage_site), (0, 1.25), c='darkgreen', linestyle='--', label='CS')


    #adjust lim to fit pos labels
    plt.ylim((-0.15,1.05))
    plt.xlim((0, seq_length+1))

    plt.ylabel('Probability')
    plt.xlabel('Protein sequence')
    if title is not None:
        plt.title(title)

    plt.legend(loc='upper right')
    plt.tight_layout()

    return fig



def make_sequence_txt_output(identifier:str, aa_sequence:str, viterbi_path: np.ndarray, marginal_probs:np.ndarray, kingdom:str, file_path:str):
    '''Write the per-position properties (viterbi label + marginal probs) to txt file
    1 line - 1 position
    Assume we are writing preprocessed viterbi paths + marginal probs,
    so we only have to write a select set of indices:
    LABEL_IDS_TO_PLOT and 0

    '''
    label_ids_to_write = [0] + LABEL_IDS_TO_PLOT
    with open(file_path, 'w') as f:
        f.write(f'# Name={identifier}\n')

        #make the column headers
        prob_col_headers = ['Other'] + [IDX_LABEL_MAP[i] for i in LABEL_IDS_TO_PLOT] # override label of pos 0
        f.write('# pos\taa\tlabel\t'+ '\t'.join(prob_col_headers)+'\n')

        for i in range(len(aa_sequence[:70])): #need to iterate over truncated aa_sequence to not exceed indexing: viterbi paths are always padded to 70
            # first 3 cols
            f.write('\t'.join([
                                str(i+1),
                                aa_sequence[i], 
                                IDS_TO_PLOT_LABEL[viterbi_path[i]]
                                ]
                                )
                                +'\t')
            # marginal probs
            probs_at_pos = marginal_probs[i]
            probs = [f'{p:.6f}' for p in probs_at_pos[label_ids_to_write]]
            f.write('\t'.join(probs)+'\n')



def make_plots_one_sample(marginal_probs, viterbi_path, aa_sequence, cleavage_site, identifier, kingdom, txt_path, png_path='', eps_path=''):
    '''Creates and saves the single-sequence output files for one sample'''

    if txt_path != '':
        make_sequence_txt_output(identifier, aa_sequence,viterbi_path,marginal_probs, kingdom, txt_path)

    if png_path != '' and eps_path != '':
        fig = sequence_plot(marginal_probs, viterbi_path, aa_sequence, cleavage_site, title='SignalP 6.0 prediction: ' +identifier)
        fig.savefig(png_path, format='png')
        fig.savefig(eps_path, format='eps')
    elif png_path != '':
        fig = sequence_plot(marginal_probs, viterbi_path, aa_sequence, cleavage_site, title='SignalP 6.0 prediction: ' +identifier)
        fig.savefig(png_path, format='png')
    elif eps_path !='':
        fig = sequence_plot(marginal_probs, viterbi_path, aa_sequence, cleavage_site, title='SignalP 6.0 prediction: ' +identifier)
        fig.savefig(eps_path, format='eps')

    plt.clf()


def make_plots_multiprocess(identifiers: List[str], sequences: List[str], global_probs: np.ndarray, cleavage_sites: np.ndarray, 
                            marginal_probs: np.ndarray, viterbi_paths:np.ndarray, output_dir:str, out_format='txt', kingdom='Other', write_procs = 8):
    '''
    Creates the single-sequence output files. Uses multiprocessing
    to speed up matplotlib plotting (bottleneck for large jobs)
    txt_only skips plot generation, further speeding up output.
    '''
    # set up the iterables for the save paths
    txt_paths = []
    eps_paths = []
    png_paths = []



    for idx, identifier in enumerate(identifiers):
        
        #prepare the file names.
        masked_identifier = re.sub('[^A-Za-z0-9]','_',identifiers[idx])

        # Always make txt.
        txt_fname = 'output_%s_plot.txt' % (masked_identifier)
        txt_path = os.path.join(output_dir, txt_fname)

        if out_format == 'none':
            txt_path = ''
            eps_path = ''
            png_path = ''

        elif out_format =='txt':
            eps_path = ''
            png_path = ''

        elif out_format == 'png':
            png_fname = 'output_%s_plot.png' % (masked_identifier)
            png_path = os.path.join(output_dir, png_fname)
            eps_path = ''

        elif out_format == 'eps':
            eps_fname = 'output_%s_plot.eps' % (masked_identifier)
            eps_path = os.path.join(output_dir, eps_fname)
            png_path = ''

        else:
            txt_fname = 'output_%s_plot.txt' % (masked_identifier)
            eps_fname = 'output_%s_plot.eps' % (masked_identifier)
            png_fname = 'output_%s_plot.png' % (masked_identifier)

            txt_path = os.path.join(output_dir, txt_fname)
            eps_path = os.path.join(output_dir, eps_fname)
            png_path = os.path.join(output_dir, png_fname)

        txt_paths.append(txt_path)
        eps_paths.append(eps_path)
        png_paths.append(png_path)


    kingdoms = [kingdom] * len(identifiers)
    
    # if we don't write output files, just skip.
    if out_format != 'none':

        fn_args = tqdm(list(zip(marginal_probs, viterbi_paths, sequences, cleavage_sites, identifiers, kingdoms, txt_paths, png_paths, eps_paths)), total=len(kingdoms), desc='Writing files')
        if write_procs >1:
            # set up multiprocessing
            pool = multiprocessing.Pool(min(write_procs, os.cpu_count()))
            pool.starmap(make_plots_one_sample, fn_args, chunksize=1)
        else:
            print('signle proc.')
            for args in fn_args:
                make_plots_one_sample(*args)
        

    return eps_paths, png_paths, txt_paths



