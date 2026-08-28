'''
Utility functions for SignalP 6.0
Input parsing, tokenization and output postprocessing
'''
from itertools import groupby
from collections import defaultdict
import numpy as np
import torch

GLOBAL_LABEL_DICT = {0:'NO_SP', 1:'SP', 2:'LIPO', 3:'TAT', 4:'TATLIPO', 5:'PILIN'}
TOKENIZER_VOCAB = {
                    '[PAD]': 0,
                    '[UNK]': 1,
                    '[CLS]': 2,
                    '[SEP]': 3,
                    '[MASK]':4,
                    'L':     5,
                    'A':     6,
                    'G':     7,
                    'V':     8,
                    'E':     9,
                    'S':    10,
                    'I':    11,
                    'K':    12,
                    'R':    13,
                    'D':    14,
                    'T':    15,
                    'P':    16,
                    'N':    17,
                    'Q':    18,
                    'F':    19,
                    'Y':    20,
                    'M':    21,
                    'H':    22,
                    'C':    23,
                    'W':    24,
                    'X':    25,
                    'U':    26,
                    'B':    27,
                    'Z':    28,
                    'O':    29,
                    # The choice of this has virtually no effect on the prediction output, see randomization experiment.
                    'eukarya': 30,
                    'other': 30, # if other is chosen, any is fine. just use euk.
                    'archaea': 31,
                    'positive':32,
                    'negative':33,
                    }

#handle non-vocab tokens: all X
TOKENIZER_VOCAB = defaultdict(lambda: 25, TOKENIZER_VOCAB)


def tokenize_sequence(amino_acids: str, kingdom: str):
    '''
    Hard-coded version of the HF tokenizer.
    Possible because AA vocabs are small, removes dependency on transformers.
    '''
    tokenized =  [TOKENIZER_VOCAB[x.upper()] for x in amino_acids]
    tokenized = [TOKENIZER_VOCAB['[CLS]'], TOKENIZER_VOCAB[kingdom]] + tokenized + [TOKENIZER_VOCAB['[SEP]']]

    return tokenized

def parse_fasta(fastafile: str):
    '''
    Parses fasta file into lists of identifiers and sequences.
    Can handle multi-line sequences and empty lines.
    Appends numeration if there are duplicate identifiers.   
    '''
    ids = []
    seqs = []
    with open(fastafile, 'r') as f:

        id_seq_groups = (group for group in groupby(f, lambda line: line.startswith(">")))

        for is_id, id_iter in id_seq_groups:
            if is_id: # Only needed to find first id line, always True thereafter
                ids.append(next(id_iter).strip())
                seqs.append("".join(seq.strip() for seq in next(id_seq_groups)[1]))
        
    #truncate and un-duplicate identifiers
    #ids = [x[:80] for x in ids] #no more truncation.
    already_seen = defaultdict(int)
    outnames = []
    for i in ids:
        #replace special characters        
        #i = i.replace(' ','_') 
        already_seen[i] += 1
        if already_seen[i]>1:
            outnames.append(i + '_' + str(already_seen[i]))
        else:
            outnames.append(i)
    ids = outnames
    #remove whitespace in fasta
    seqs = [x.replace(' ', '') for x in seqs]

    return ids, seqs


# Eukarya summing indices
all_n_region = [3,9,16,17,23,24,31]
all_h_region = [4,10,18,25]
all_c_region = [5,11,19,26] # Sec/SPI c, Sec/SPII c, Tat/SPI c, Tat/SPII c
all_mature = [12,27,32,33] #not really all mature. This is just additional for Euk
# Other summing indices
LABEL_IDS_OTHER = [0,1,2,6,7,8,13,14,15,20,21,22,28,29,30,34,35,36]

def postprocess_probabilities(all_probs: np.ndarray, kingdom_id: str):
    '''
    Process the marginal probabilities to get cleaner outputs.
    Eukarya: sum all other classes into Sec/SPI
    Sec/SPII: sum h region and lipobox
    Sum all I,M,O probabilities into 0 (for plotting)
    '''
    all_probs_out = all_probs.copy()

    # set 0 in out, then fill in again.
    if kingdom_id  == 'eukarya':
        all_probs_out[:,:,all_n_region] = 0
        all_probs_out[:,:,all_h_region] = 0
        all_probs_out[:,:,all_c_region] = 0
        all_probs_out[:,:,all_mature] = 0

        all_probs_out[:,:,3] = all_probs[:,:,all_n_region].sum(axis=-1)
        all_probs_out[:,:,4] = all_probs[:,:,all_h_region].sum(axis=-1)
        all_probs_out[:,:,5] = all_probs[:,:,all_c_region].sum(axis=-1)
        all_probs_out[:,:,6] = all_probs[:,:,all_mature].sum(axis=-1)

    else:
        # no c region for SPII
        h_sec_probs= all_probs[:,:,[10,11]].sum(axis=-1) # merge Sec/SPII h + c probs
        h_tat_probs= all_probs[:,:,[25,26]].sum(axis=-1) # merge Tat/SPII h + c probs

        all_probs_out[:,:,11] = 0 
        all_probs_out[:,:,26] = 0
        all_probs_out[:,:,10] = h_sec_probs
        all_probs_out[:,:,25] = h_tat_probs

    # merge all OTHER probs into 0
    all_probs_out[:,:,0] = all_probs[:,:, LABEL_IDS_OTHER].sum(axis=-1)
    
    return all_probs_out

def postprocess_viterbi(viterbi_paths: np.ndarray, kingdom_id: str):
    '''Replace labels in viterbi paths for simplification'''
    viterbi_paths_out = viterbi_paths.copy()
    if kingdom_id == 'eukarya':
        #replace all n with Sec/SPI n
        mask = np.isin(viterbi_paths, all_n_region)
        viterbi_paths_out[mask] =3
        #replace all h with Sec/SPI h
        mask = np.isin(viterbi_paths, all_h_region)
        viterbi_paths_out[mask] =4
        #replace all c with Sec/SPI c
        mask = np.isin(viterbi_paths, all_c_region)
        viterbi_paths_out[mask] =5
        # replace modified cysteine with OTHER
        mask = np.isin(viterbi_paths, all_mature)
        viterbi_paths_out[mask] = 0  
    else:
        #lipobox masking - replace lipobox with h
        viterbi_paths_out[viterbi_paths==11] = 10
        viterbi_paths_out[viterbi_paths==26] = 25

    # merge I,M,O states into OTHER
    viterbi_paths_out[np.isin(viterbi_paths, LABEL_IDS_OTHER)] = 0 

    return viterbi_paths_out

            


def model_inputs_from_fasta(fastafile: str, kingdom_id: str):
    '''Parse a fasta file to input id + mask tensors.
    Pad all seqs to full length (73: 70+special tokens).
    traced model requires that, it was traced with 73 as fixed length'''


    identifiers, sequences = parse_fasta(fastafile)

    if len(identifiers) == 0:
        return [], [], torch.LongTensor(), torch.LongTensor()

    #truncate
    input_ids = [x[:70] for x in sequences]
    input_ids =  [tokenize_sequence(x, kingdom_id) for x in input_ids]
    input_ids = [x + [0] * (73-len(x)) for x in input_ids]

    input_ids = np.vstack(input_ids)
    input_mask = (input_ids>0) *1    

    identifiers = [x.lstrip('>') for x in identifiers]

    return identifiers, sequences, torch.LongTensor(input_ids), torch.LongTensor(input_mask)




def get_cleavage_sites(tagged_seqs: np.ndarray, cs_tokens = [5, 11, 19, 26, 31]):
    '''Convert sequences of tokens to the indices of the cleavage site.
    Inputs:
        tagged_seqs: (batch_size, seq_len) integer array of position-wise labels , with "C" tags for cleavage sites
        cs_tokens = label tokens that indicate a cleavage site
    Returns:
        cs_sites: (batch_size) integer array of position that is a CS. -1 if no SP present in sequence.
    '''
    def get_last_sp_idx(x: np.ndarray) -> int:
        '''Func1d to get the last index that is tagged as SP. use with np.apply_along_axis. '''
        sp_idx = np.where(np.isin(x, cs_tokens))[0]
        if len(sp_idx)>0:
            max_idx = sp_idx.max()+1
        else:
            max_idx = -1
        return max_idx

    cs_sites = np.apply_along_axis(get_last_sp_idx, 1, tagged_seqs)
    return cs_sites

def get_cleavage_site_probs(cleavage_sites, viterbi_paths, marginal_probs):
    '''Get the cleavage site probability (= probability of preceding token that caused the CS annotation)
    get the label idx of the found CS from the viterbi path,
    use it to index into marginal probs
    '''
    #need to subtract 1 from cs again for indexing
    last_sp_pos_label = viterbi_paths[np.arange(viterbi_paths.shape[0]), cleavage_sites-1]
    probs =  marginal_probs[np.arange(viterbi_paths.shape[0]),cleavage_sites-1, last_sp_pos_label]
    #mask out 
    probs[cleavage_sites==-1] = 0
    return probs


# States for types.
# tuple (SP_STATES, OTHER_STATES)
TYPE_STATES = {
    0: ([],[]), # NO_SP , don't need any states for path imputation
    1: ([3,4,5], [6,7,8]), # SP
    2: ([9,10,11,12], [13,14,15]), # LIPO
    3: ([16,17,18,19], [20,21,22]), # TAT
    4: ([23,24,25, 26, 27], [28,29,30]), # TATLIPO
    5: ([31, 32, 33], [34,35,36]), # PILIN
}

def check_viterbi_type(viterbi_path: np.ndarray) -> int:
    '''Returns the global type of the viterbi path.'''
    first_label =  viterbi_path[0]
    if first_label in [0,1,2]:
        return 0

    for type_label, state_labels in TYPE_STATES.items():
        sp_states = state_labels[0]
        if first_label in sp_states:
            return type_label
        
    raise NotImplementedError(f'A viterbi path of a SP cannot start with {first_label}. Model must have failed internally, report a bug to developer.')

def fix_twin_arginine(viterbi_path: np.ndarray, marginal_probabilities: np.ndarray, pred_type:int) -> np.ndarray:
    '''
    Ensures that a twin-RR motif of two adjacent R states is predicted, followed by exactly one N.
    '''
    viterbi_path = viterbi_path.copy()

    arginine_token = 17 if pred_type == 3 else 24
    n_token = 16 if pred_type == 3 else 23
    h_token = 18 if pred_type == 3 else 25
    
          
    # Step 1: Ensure there is a twin-RR motif predicted.
    twin_arg_positions = np.where(viterbi_path==arginine_token)[0]
    if len(twin_arg_positions) == 1:
        # find the neighboring pos with highest R prob and relabel to R
        # This code will fail in the (unlikely and broken) case that R is the last position. Rather keep bad prediction than crash.
        try:
            prob_before =  marginal_probabilities[twin_arg_positions[0]+1, arginine_token]
            prob_after =  marginal_probabilities[twin_arg_positions[0]-1,arginine_token]
            if prob_before>prob_after:
                viterbi_path[twin_arg_positions[0]+1] = arginine_token
            else:
                viterbi_path[twin_arg_positions[0]-1] = arginine_token
        except:
            pass

    elif len(twin_arg_positions) > 2:
        # find the most probable R.
        prob_order = np.argsort(marginal_probabilities[twin_arg_positions, arginine_token])[::-1]
        arginines_ordered = twin_arg_positions[prob_order]
        best_arginine =  arginines_ordered[0]
        # find the best one that is adjacent to the best.
        is_adjacent = np.abs(arginines_ordered[1:] - best_arginine) ==1
        if is_adjacent.sum() == 0:
            # this should not happen - no R predicted next to best R. Just use a really basic default rule.
            second_best_arginine = best_arginine - 1
        else:
            second_best_arginine = arginines_ordered[1:][is_adjacent][0]

        # all remaining Rs are set to N.
        viterbi_path[twin_arg_positions] = n_token
        viterbi_path[best_arginine] = arginine_token
        viterbi_path[second_best_arginine] = arginine_token


    # Step 2: ensure exactly one N comes after the RR. Replaces excess with H.
    last_valid_n_position = np.where(viterbi_path==arginine_token)[0][1] + 1
    viterbi_path[last_valid_n_position] = n_token


    n_pos = np.where(viterbi_path==n_token)[0]
    n_pos = n_pos[n_pos>last_valid_n_position]
    viterbi_path[n_pos] = h_token
     
    return viterbi_path

def resolve_viterbi_marginal_conflicts(global_probs: np.ndarray, marginal_probs: np.ndarray, cleavage_sites: np.ndarray, viterbi_paths: np.ndarray) -> None:
    '''
    This function needs to be called after all other pre-processings steps.

    Viterbi decoding and marginal probablities sometimes can show conflicting predictions:
    - marginal probabilities for SP too low, predicted as other, but viterbi path with CS found  
    - both predict a SP, but the type is different.

    The marginal prediction takes precedence over the viterbi path. Conflicts are resolved by
    - Marginal probs too low: `viterbi_paths` is masked to Other, no CS reported.  
    - Different type: Replace the viterbi path with the argmax of the marginal probabilities.

    We take advantage of the already predicted SP type to mask out the marginals
    of all non-associated states when computing the argmax. In case the argmax path does not follow
    the predefined grammar (only one R, multiple cysteines, mixed regions), we resolve these issues 
    heuristically.

    The whole operation is wrapped in a try-except block, as reporting broken predictions is more useful
    than the prediction jobs crashing completely.

    Modifies input arrays `viterbi_paths` and `cleavage_sites` inplace.
    '''
    failed_resolve = []
    for i in range(global_probs.shape[0]):

        pred_type = np.argmax(global_probs[i])
        viterbi_type = check_viterbi_type(viterbi_paths[i]) # check in which CRF branch the viterbi path is.

        try:
            # Case 1: Marginal probabilities too low.
            # global label NO_SP and cleavage site predicted.
            if (pred_type == 0) & (cleavage_sites[i] != -1):
                viterbi_paths[i, :] = 0 # set viterbi path to NO_SP
                cleavage_sites[i] = -1 # don't predict a CS.

            
            # Case 2: Viterbi path has different type than global prediction.        
            elif pred_type != viterbi_type:
                
                # 1. Subset probability channels to all states associated with this type.
                sp_states, other_states = TYPE_STATES[pred_type]
                type_marginal_probs = np.zeros_like(marginal_probs[i])
                type_marginal_probs[:,sp_states] = marginal_probs[i][:, sp_states]
                type_marginal_probs[:,0] = marginal_probs[i][:,other_states].sum(axis=-1)


                # 2. Get the argmax path and substitute the viterbi path with it.
                marginal_region_preds = np.argmax(type_marginal_probs, -1)
                sp_idx = np.where(np.isin(marginal_region_preds, [5, 10, 19, 25, 31]))[0]
                viterbi_paths[i] = marginal_region_preds

                # 3. Get the cleavage site.
                if pred_type == 2 or pred_type == 4:
                    # SPII case: because the model predicts the modified cysteine at CS+1, we can
                    # find the position with the highest probability and call the CS there.
                    modified_cys_probs = marginal_probs[i,:, [12,27]].sum(axis =0)
                    putative_modified_cys_pos = np.argmax(modified_cys_probs) # find the position with maximum C+1 probability.
                    cleavage_sites[i] = putative_modified_cys_pos # no need for +1, as is already.
                    # also need to fix the viterbi path.
                    viterbi_paths[i, putative_modified_cys_pos] = 12 if pred_type == 2 else 27

                else:
                    # Default case: set CS to last sp position in surrogate viterbi path.
                    cleavage_sites[i] = sp_idx.max() +1

                # 4. Fix potential issues in the viterbi path.
                
                # Everything past the CS needs to be other.
                if pred_type == 2  or pred_type == 4:
                    viterbi_paths[i, cleavage_sites[i]+1:] = 0
                else:
                    viterbi_paths[i, cleavage_sites[i]:] = 0

                # At region transition positions, OTHER can be maximum.
                # e.g. n=0.2, h=0.35, O=0.45 --> argmax is O although it's below 50%, but we want h.
                # do argmax with OTHER probs masked out up to CS.
                type_marginal_probs[:cleavage_sites[i], 0] = 0
                viterbi_paths[i, :cleavage_sites[i]] = np.argmax(type_marginal_probs[:cleavage_sites[i]], -1)


                # for TAT SPs, argmax might fail to give good twin arginines.
                if pred_type == 3 or pred_type == 4:
                    fixed_path = fix_twin_arginine(viterbi_paths[i], marginal_probs[i], pred_type)
                    viterbi_paths[i] = fixed_path
                    
                # for n, h, c, sometimes the argmax path can be e.g.
                # NNNNNNNNNHHHHHHHHHNNHHHHHHHHCCHCCCCOOOOOOO
                # Use a heuristic fixing order to decide:

                # After the first h, there may be no more N
                # After the first c, there may be no more H

                # before the last N, everything is N.


                last_n = np.where(np.isin(viterbi_paths[i], all_n_region))[0].max()
                n_label = viterbi_paths[i, last_n] # retrieve the correct n label for this type.
                twin_arg_label = 17 if pred_type == 3 else 24
                twin_args = np.where(viterbi_paths[i] == twin_arg_label)[0]
                viterbi_paths[i, :last_n+1] =  n_label
                viterbi_paths[i, twin_args] =  twin_arg_label

                # before the last H and after the last N, everything is H.
                last_h = np.where(np.isin(viterbi_paths[i], all_h_region))[0].max()
                h_label = viterbi_paths[i, last_h]
                viterbi_paths[i, last_n+1:last_h+1] = h_label

        except:
            failed_resolve.append(i)

    return failed_resolve



            
            