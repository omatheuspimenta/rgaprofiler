'''
Main prediction script.
Used by webserver, and can also be used from CLI.
Due to reliance on TorchScript, this can be effectively run witout
any dependencies - all we need is the compiled model.

The script works on a fasta file.

Args:
    --fastafile 
    --output_path
    --format
    --organism

Example :
    python3 predict.py --fastafile path/to/file.fasta --output output/dir/ --organism Archaea


results are collected in results_summary dict to be dumped as json.
 -general job info in "INFO"
 -individual outputs in "SEQUENCES


'''
# torch gives deprecation warnings. Suppress for now, once the deprecation is here freeze the torch version.
import warnings
warnings.filterwarnings("ignore", category=UserWarning)

import torch
import json
import argparse
import os
import numpy as np
from datetime import datetime
from shutil import make_archive, move
from tqdm import tqdm

from .sequential_model import predict_sequential
from .make_sequence_plot import make_plots_multiprocess
from .make_output_files import make_output_table, make_gff_file, make_region_gff_file, make_fasta
from .utils import model_inputs_from_fasta, get_cleavage_sites, postprocess_probabilities, postprocess_viterbi, get_cleavage_site_probs, resolve_viterbi_marginal_conflicts
import traceback
from . import __version__

BASE_DIR = os.path.dirname(__file__)
MODEL_DIR = os.path.join(BASE_DIR, 'model_weights')

def predict(model: torch.nn.Module, input_ids: torch.Tensor, input_mask: torch.Tensor, batch_size: int=10):
    '''Cut batches from tensors and process batch-wise'''
    #process as batches
    all_global_probs = []
    all_marginal_probs = []
    all_viterbi_paths = []
    b_start = 0
    b_end = batch_size

    device = model.crf.transitions.device #Get the device of the jitted model. Any weight works.
    model.to(device)

    pbar = tqdm(total=len(input_ids), desc='Predicting', unit='sequences')

    while b_start < len(input_ids):
        with torch.no_grad():
            ids = input_ids[b_start:b_end,:].to(device)
            mask =  input_mask[b_start:b_end,:].to(device)
            global_probs, marginal_probs, viterbi_paths = model(ids, mask)
            
            #this did not work in tracing yet.
            viterbi_paths = torch.nn.functional.pad(viterbi_paths,(0, 70-viterbi_paths.shape[1]), value=-1)

            all_global_probs.append(global_probs.cpu().numpy())
            all_marginal_probs.append(marginal_probs.cpu().numpy())
            all_viterbi_paths.append(viterbi_paths.cpu().numpy())

        b_start = b_start + batch_size
        b_end = b_end + batch_size
        pbar.update(ids.shape[0])

    pbar.close()
    all_global_probs = np.concatenate(all_global_probs)
    all_marginal_probs = np.concatenate(all_marginal_probs)
    all_viterbi_paths = np.concatenate(all_viterbi_paths)

    return all_global_probs, all_marginal_probs, all_viterbi_paths



def add_sequences_to_summary(summary:dict, identifiers:list, global_probs:np.ndarray, cleavage_sites:np.ndarray, 
                                                    cs_probs: np.ndarray, kingdom:str, eps_save_paths, png_save_paths, txt_save_paths):
    '''Adds each predicted sequence to the JSON output. Paths need to be converted before adding to the JSON.'''

    if kingdom == 'eukarya':
        p_no =  global_probs[:,0, None]
        p_sp = global_probs[:,1:, None].sum(axis=1)
        global_probs = np.hstack([p_no, p_sp])
        prednames = ['Other','Signal Peptide (Sec/SPI)']
    else:
        prednames = ['Other', 'Signal Peptide (Sec/SPI)', 'Lipoprotein signal peptide (Sec/SPII)', 'TAT signal peptide (Tat/SPI)',
                     'TAT Lipoprotein signal peptide (Tat/SPII)', 'Pilin-like signal peptide (Sec/SPIII)']

    pred_label_id = np.argmax(global_probs, axis=1)

    for idx, identifier in enumerate(identifiers):
        seq_dict = {
            'Name': identifier,
            'Plot_eps': eps_save_paths[idx],
            'Plot_png': png_save_paths[idx],
            'Plot_txt': txt_save_paths[idx],
            'Prediction': prednames[pred_label_id[idx]],
            'Likelihood': [round(x.item(), 4) for x in global_probs[idx]],
            'Protein_types': prednames

        }
        if cleavage_sites[idx] == -1:
            seq_dict['CS_pos'] = ''
        else:
            seq_dict['CS_pos'] = f"Cleavage site between pos. {cleavage_sites[idx]} and {cleavage_sites[idx]+1}. Probability {cs_probs[idx].item():.6f}"

        summary['SEQUENCES'][identifier] = seq_dict

    return summary



def main():
    parser = argparse.ArgumentParser('SignalP 6.0 Signal peptide prediction tool', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--fastafile', '-ff' ,'-fasta', type=str, help='Amino acid sequences to predict in FASTA format.', required=True)
    parser.add_argument('--output_dir', '-od', type=str, help='Path at which to save the output files. Will be created if not existing already.', required=True)
    parser.add_argument('--format', '-fmt', '-f', '-format', type=str, default='txt', const='txt', nargs='?', choices=['txt', 'png', 'eps', 'all', 'none'],
                        help='Type of single-sequence output files to be created. `txt` produces tabular output, `png`, `eps` and `all` additionally produce figures. `none` only creates prediction summary files.')
    parser.add_argument('--organism', '-org', type=str, help='The organism group of origin of the sequences. `eukarya`, `euk` limit predictions to Sec/SPI.', default='other', const='other', nargs='?', choices=['eukarya', 'other',  'euk'])
    parser.add_argument('--mode', '-m', default='fast', const='fast', nargs='?', choices=['fast', 'slow', 'slow-sequential'], help='Which prediction model to use. Modes might have to be installed manually.')
    parser.add_argument('--bsize', '-bs', '-batch', type=int, default=10, help='Batch size for prediction. Use to adjust performance to your system memory.')
    parser.add_argument('--write_procs', '-wp', type=int, default=8, help='Number of parallel processes used for writing ouput files. Use to adjust performance to your system memory.')
    parser.add_argument('--torch_num_threads', '-tt', type=int, default=8, help='Number of threads used by PyTorch for neural network execution.')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    # Advanced.
    parser.add_argument('--skip_resolve', action='store_true', help="Skip resolving conflicts of the model's Viterbi path and marginal probabilities. In general, should not be activated as it might cause inconsistent prediction behaviour.")
    parser.add_argument('--model_dir','-md', type=str, help='Path to directory that contains the model files. Use only when not following default installation instructions.', default=MODEL_DIR)
    

    args = parser.parse_args()

    torch.set_num_threads(args.torch_num_threads)


    slow_model_path = os.path.join(args.model_dir,'ensemble_model_signalp6.pt')
    fast_model_path = os.path.join(args.model_dir, 'distilled_model_signalp6.pt')
    sequential_model_path = os.path.join(args.model_dir, 'sequential_models_signalp6')

    # For backwards compatibility with SP5.
    if args.organism == 'euk':
        setattr(args, 'organism', 'eukarya')

    if args.fastafile is None:
        raise ValueError('Argument "--fastafile" needs to be specified.')
    if args.output_dir is None:
        raise ValueError('Argument "--output_dir" needs to be specified.')
    if not os.path.isfile(args.fastafile):
        raise ValueError(f'Specified file path {args.fastafile} is not a existing file.')
    if not os.path.exists(args.output_dir):
        try:
            os.makedirs(args.output_dir)
        except:
            raise ValueError(f'{args.output_dir} does not exist and cannot be created.')

    # Parse inputs
    identifiers, aa_sequences, input_ids, input_mask = model_inputs_from_fasta(args.fastafile, args.organism)

    #Prepare output JSON object to be filled
    result_summary = {
        "INFO": {
            "failedjobs": 0,
            "size": len(identifiers)
        },
        "CSV_FILE": "",
        "MATURE_FILE": "",
        "GFF_FILE": "",
        "ZIP_FILE": "",
        "SEQUENCES": {},
        "FORMAT": args.format,
        "ORG": args.organism,
        "MODE": args.mode,
    }


    # TODO add raise statements instead of writing error.json for local version.
    # this only made sense as the server backend.           
    # validate inputs, break from main if not parseable
    if len(identifiers) != len(aa_sequences):
        with open(os.path.join(args.output_dir,'errors.json'), 'w+') as outf:
            outf.write("{ parse_error: %r}" % ("Could not parse input!" +
                        " Please report an error.") )
        result_summary['INFO']['failedjobs'] = 1
        raise ValueError(f'Fasta file {args.fastafile} could not be parsed correctly. Different number of headers and sequences.')
        return None


    # Check the parser output whether there's something wrong
    # If the parser gives us no data, write error file to output dir
    if len(identifiers) < 1:
        with open(os.path.join(args.output_dir,'errors.json'), 'w+') as outf:
            outf.write("{ parse_error: %r}" % ("Could not parse input!" +
                        " Check that data file is not empty / is in fasta format.") )
        result_summary['INFO']['failedjobs'] = 1
        raise ValueError(f'Fasta file {args.fastafile} could not be parsed correctly. Check that file is not empty / is in fasta format.')
        return None




    # load model and run data

    if args.mode =='slow':
        try:
            model = torch.jit.load(slow_model_path)
        except RuntimeError:
            raise RuntimeError('The model checkpoint is configured to run on GPU, but no CUDA device was found on the current system. Revert back to CPU or fix NVIDIA drivers.')
        except ValueError:
            # check if slow-sequential is available, else error.
            if os.path.exists(sequential_model_path):
                setattr(args, 'mode', 'slow-sequential')
                # print('Mode `slow` not found, running `slow-sequential` instead.')
            else:
                raise FileNotFoundError(f"Slow mode requires the full model to be installed at {slow_model_path}. Missing from this installation or incorrectly configured.")
    elif args.mode == 'slow-sequential':
        if not os.path.exists(sequential_model_path):
            raise FileNotFoundError(f'Sequential mode requires the models to be installed at {sequential_model_path}. Missing from this installation.')
    else:
        try:
            model = torch.jit.load(fast_model_path)
        except RuntimeError:
            raise RuntimeError('The model checkpoint is configured to run on GPU, but no CUDA device was found on the current system. Revert back to CPU or fix NVIDIA drivers.')
        except ValueError:
            raise FileNotFoundError(f"Fast mode requires model to be installed at {fast_model_path}. Missing from this installation.")

    if args.mode == 'slow-sequential':
        global_probs, marginal_probs, viterbi_paths =  predict_sequential(sequential_model_path, input_ids, input_mask, batch_size=args.bsize)
    else: 
        global_probs, marginal_probs, viterbi_paths =  predict(model, input_ids, input_mask, batch_size=args.bsize)
    

    ## Post-process the raw model outputs to get all predictions in the desired format.

    # 1. Infer the cleavage sites from the viterbi paths.
    cleavage_sites =  get_cleavage_sites(viterbi_paths) #this uses labels before post-processing.

    # 2. Simplify the marginal probabilities and merge classes into Sec/SPI for eukarya.
    marginal_probs = postprocess_probabilities(marginal_probs, args.organism)

    # 3. Same simplifications for viterbi paths.
    viterbi_paths = postprocess_viterbi(viterbi_paths, args.organism)

    # 4. Merge global class probabilities into Sec/SPI for eukarya.
    if args.organism == 'eukarya':
        global_probs[:,1] = global_probs[:,1:].sum(axis=1)
        global_probs[:,2:] = 0

    # 5. Resolve edge case discrepancies between viterbi decoding and marginal probabilities.
    if not args.skip_resolve:
        failed_resolve = resolve_viterbi_marginal_conflicts(global_probs, marginal_probs, cleavage_sites, viterbi_paths)
        for i in failed_resolve:
            print(f'Unknown behaviour encountered for sequence {i+1}. Please check outputs: {identifiers[i]}')

    # 6. Retrieve the probability of each predicted cleavage site (=probability of preceding label)
    cs_probs = get_cleavage_site_probs(cleavage_sites, viterbi_paths, marginal_probs)

    # 7. clip probabilities to 0,1 before writing outputs
    marginal_probs = np.clip(marginal_probs, 0, 1)
    global_probs = np.clip(global_probs, 0, 1)
    cs_probs = np.clip(cs_probs, 0, 1)


    # write job output files
    tsv_output_file  = os.path.join(args.output_dir, 'prediction_results.txt')
    make_output_table(identifiers, global_probs, cleavage_sites, cs_probs, args.organism, tsv_output_file)
    result_summary['CSV_FILE'] = tsv_output_file

    gff_output_file  = os.path.join(args.output_dir, 'output.gff3')
    make_gff_file(identifiers, global_probs, cleavage_sites, args.organism, gff_output_file)
    result_summary['GFF_FILE'] = gff_output_file

    gff_output_file  = os.path.join(args.output_dir, 'region_output.gff3')
    make_region_gff_file(identifiers, global_probs, viterbi_paths, cleavage_sites, args.organism, gff_output_file)
    result_summary['REGION_GFF_FILE'] = gff_output_file

    fasta_output_file  = os.path.join(args.output_dir, 'processed_entries.fasta')
    make_fasta(identifiers, aa_sequences, fasta_output_file, cleavage_sites)
    result_summary['MATURE_FILE'] = fasta_output_file



    # write single-sequence output files
    eps_save_paths, png_save_paths, txt_save_paths = make_plots_multiprocess(identifiers, aa_sequences, global_probs, cleavage_sites, marginal_probs, viterbi_paths, args.output_dir, out_format=args.format, kingdom=args.organism, write_procs=args.write_procs)
    add_sequences_to_summary(result_summary, identifiers, global_probs, cleavage_sites, cs_probs, args.organism, eps_save_paths, png_save_paths, txt_save_paths)



    # Make JSON output
    with open(os.path.join(args.output_dir,'output.json'),'w+') as outf:
        json.dump(result_summary,outf,indent=2)

     
    
if __name__ == "__main__":
    main()
    
