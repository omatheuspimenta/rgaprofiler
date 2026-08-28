'''
Utility functions to generate all tabular output files.
'''
import numpy as np
from typing import List, Tuple
from datetime import datetime

GLOBAL_LABEL_DICT = {0:'OTHER', 1:'SP', 2:'LIPO', 3:'TAT', 4:'TATLIPO', 5:'PILIN'}

def make_output_table(identifiers, global_probs, cleavage_sites, cleavage_sites_probs, kingdom, output_file_path):
    '''Make the .tsv tabular output for all sequences'''
    timestamp = datetime.now().strftime('%Y%m%d%H%M%S')

    # get the global preds, post-processing for Eukarya
    if kingdom == 'eukarya':
        p_no =  global_probs[:,0, None] # keep a dummy dim to put array back together
        p_sp = global_probs[:,1:, None].sum(axis=1)
        global_probs = np.hstack([p_no, p_sp])
        colnames = ['ID', 'Prediction', 'OTHER', 'SP(Sec/SPI)', 'CS Position']
    else:
        colnames = ['ID', 'Prediction', 'OTHER', 'SP(Sec/SPI)', 'LIPO(Sec/SPII)', 'TAT(Tat/SPI)', 'TATLIPO(Tat/SPII)', 'PILIN(Sec/SPIII)', 'CS Position']
    
    pred_label_id = np.argmax(global_probs, axis=1)

    with open(output_file_path, 'w') as f:
        #write the headers
        f.write(f'# SignalP-6.0\tOrganism: {kingdom.lower().capitalize()}\tTimestamp: {timestamp}\n')
        f.write('# '+'\t'.join(colnames) +'\n')

        for idx, identifier in enumerate(identifiers):
            #format everything for writing
            prediction = GLOBAL_LABEL_DICT[pred_label_id[idx]]
            probs = global_probs[idx]
            probs = [f'{p:.6f}' for p in probs]
            cs = cleavage_sites[idx] #if cleavage_sites[idx] != -1 else ''
            cs_prob = cleavage_sites_probs[idx] #if cleavage_sites[idx] != -1 else ''
            cs_string = f'CS pos: {cs}-{cs+1}. Pr: {cs_prob:.4f}' if cleavage_sites[idx] !=-1 else ''
            writelist = [identifier] + [prediction] + probs + [cs_string]
            f.write('\t'.join(writelist) +'\n')




def make_gff_file(identifiers, global_probs, cleavage_sites, kingdom, output_file_path):
    '''
    Make the regular gff output, as in signalp5. Annotated features are SP types
    '''
    #sp_P14009_14KD_DAUCA	SignalP-5.0	signal_peptide	1	25	0.980592	.	.	.
    if kingdom == 'eukarya':
        p_no =  global_probs[:,0, None] # keep a dummy dim to put array back together
        p_sp = global_probs[:,1:, None].sum(axis=1)
        global_probs = np.hstack([p_no, p_sp])


    pred_label_id = np.argmax(global_probs, axis=1)

    feature_names = ['None','signal_peptide', 'lipoprotein_signal_peptide', 'signal_peptide', 'lipoprotein_signal_peptide', 'signal_peptide']

    with open(output_file_path, 'w') as f:
        f.write("## gff-version 3\n")  

        for idx, identifier in enumerate(identifiers):
            #skip no-sp sequences
            if pred_label_id[idx] != 0:
                feature_name = feature_names[pred_label_id[idx]] 
                note = 'Note=TAT' if pred_label_id[idx] in [3,4] else 'Note=Pilin' if pred_label_id[idx]==5 else '.'

                linelist = [identifier, 
                            'SignalP-6.0',
                            feature_name, 
                            '1', 
                            str(cleavage_sites[idx]), 
                            str(global_probs[idx][pred_label_id[idx]]), 
                            '.', 
                            '.',
                            note,
                            ]
                f.write('\t'.join(linelist) + '\n')



region_ids = {'n-region' : [3,9,16,23], 
              'h-region' : [4,10,18,25],
              'c-region' : [5,19],
              'lipid-modified cysteine': [12,27],
              'twin-arginine motif': [17,24]
              }


def parse_regions_from_viterbi_path(viterbi_path: np.ndarray) ->List[Tuple[str,int,int]]:
    '''Get a list of (region, start, end) tuples from a viterbi path.'''
    # get list of symbols in path
    outlist = []
    for region in region_ids.keys():
        sp_pos =  np.isin(viterbi_path, region_ids[region])
        pos_indices = np.where(sp_pos)[0]
        if len(pos_indices) >0:
            outlist.append((region, pos_indices.min()+1, pos_indices.max() +1))

    return outlist



def make_region_gff_file(identifiers, global_probs, viterbi_paths, cleavage_sites, kingdom, output_file_path):
    '''
    Make gff output for the predicted regions.
    --> multiple lines per sequence
    '''
    with open(output_file_path, 'w') as f:
        f.write("## gff-version 3\n")  

        for idx, identifier in enumerate(identifiers):
            regions = parse_regions_from_viterbi_path(viterbi_paths[idx])
            for region in regions:
                note = '.'
                linelist = [identifier,
                            'SignalP-6.0', 
                            region[0], 
                            str(region[1]), 
                            str(region[2]), 
                            '.', 
                            '.', 
                            '.',
                            note,
                            ]
                f.write('\t'.join(linelist) + '\n')


def make_fasta(identifiers: List[str], aa_sequences: List[str], output_file_path: str, cleavage_sites: List[int]):
    '''Dump SP sequences as `.fasta` . Trims signal peptide from sequence.'''
    with open(output_file_path, 'w') as f:
        for h,s, cs in zip(identifiers, aa_sequences, cleavage_sites):
            if cs != -1:

                f.write('>'+h+'\n')
                f.write(s[cs:]+'\n') # start from CS = first pos in mature seq.
