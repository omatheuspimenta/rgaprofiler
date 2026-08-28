'''
Code to run 6 SP6 models sequentially, averaging their results.
To minimize model loading times, the full data is predicted by each model,
and then averaged.
'''
import torch
from typing import List, Tuple
import os
from tqdm import tqdm

def memory_efficient_mean(tensor_list: List[torch.Tensor]) -> torch.Tensor:
    '''
    torch.stack().mean() goes OOM for large datasets,
    although the list of tensors fits memory. 
    Use this as a workaround.
    '''
    with torch.no_grad():
        sum_tensor = torch.zeros_like(tensor_list[0])
        for t in tensor_list:
            sum_tensor = sum_tensor + t

        return sum_tensor / len(tensor_list)



def get_averaged_emissions_and_probs(model_path: str, input_ids: torch.Tensor, 
                                        input_mask: torch.Tensor, 
                                        batch_size: int=100,
                                        ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    '''
    Run forward pass for each model. Average the emissions and probabilities.
    '''

    all_emissions= []
    all_global_probs = []
    all_marginal_probs = []


    for n, model in enumerate(['test_0_val_1', 'test_0_val_2', 'test_1_val_0', 'test_1_val_2', 'test_2_val_0','test_2_val_1']):
        try:
            model = torch.jit.load(os.path.join(model_path, model+'.pt'))
        except RuntimeError:
            raise RuntimeError('The model checkpoint is configured to run on GPU, but no CUDA device was found on the current system. Revert back to CPU or fix NVIDIA drivers.')
        except ValueError:
            raise FileNotFoundError('Could not find the required model. Your slow_sequential installation is probably corrupted with one or multiple files missing.')
        
        device = model.crf.transitions.device #Get the device of the jitted model. Any weight works.
        model.to(device)
        b_start = 0
        b_end = batch_size

        model_emissions = []
        model_global_probs = []
        model_marginal_probs = []

        pbar = tqdm(total=len(input_ids), desc=f'Predicting {n+1}/6',leave=True,unit='sequences')
        while b_start < len(input_ids):
            with torch.no_grad():
                ids = input_ids[b_start:b_end,:].to(device)
                mask =  input_mask[b_start:b_end,:].to(device)
                global_probs, marginal_probs, emissions = model(ids, mask)

                model_emissions.append(emissions.detach().cpu())
                model_global_probs.append(global_probs.detach().cpu())
                model_marginal_probs.append(marginal_probs.detach().cpu())


                b_start = b_start + batch_size
                b_end = b_end + batch_size
                pbar.update(ids.shape[0])
        
        pbar.close()

        # make tensors from minibatch lists
        all_emissions.append(torch.cat(model_emissions))
        all_global_probs.append(torch.cat(model_global_probs))
        all_marginal_probs.append(torch.cat(model_marginal_probs))

    # Done predicting. Average.
    mean_emissions = memory_efficient_mean(all_emissions)
    mean_global_probs = memory_efficient_mean(all_global_probs)
    mean_marginal_probs = memory_efficient_mean(all_marginal_probs)

    return mean_emissions, mean_global_probs, mean_marginal_probs


def get_viterbi_paths(model_path: torch.Tensor, emissions: torch.Tensor, input_mask: torch.Tensor):
    '''
    Load viterbi decoding module and run on provided averaged emissions.
    '''
    model = torch.jit.load(os.path.join(model_path, 'averaged_viterbi.pt'))

    device = model.transitions.device
    input_mask = input_mask.to(device)
    emissions = emissions.to(device)
    # No need to batch. Viterbi internally processes sequences one by one.
    viterbi_paths = model(emissions, input_mask)
    viterbi_paths = torch.nn.functional.pad(viterbi_paths,(0, 70-viterbi_paths.shape[1]), value=-1)

    return viterbi_paths

def make_viterbi_mask(input_mask):
    '''Trim mask and set <sep> position to 0. This is the same function we use internally in the torchscript model.'''
    input_mask = input_mask[:,2:] # remove CLS ORG --> (n,71), XXXXX<SEP>
    # shift input mask by one. the last pos will be zero, and all others are shifted to
    # the left by one. In full-length seqs, this makes the [SEP] token 0. In padded seqs.
    # this shifts the first 0 from the first [PAD] to [SEP]
    shifted_input_mask = input_mask[:,1:] # --> (n,70), XXXX<SEP>
    #add the zero vector at the end
    zeros = input_mask[:,1] * 0 #make zero vector directly from input tensor, don't compute shapes
    shifted_input_mask = torch.cat([shifted_input_mask,zeros.unsqueeze(1)], dim=1) #add dummy seq dim and concatenate along seq dim.

    #do the same to the shifted mask
    mask_out = shifted_input_mask[:,:-1]

    return mask_out

def predict_sequential(model_path: str, input_ids: torch.Tensor, input_mask: torch.Tensor, batch_size: int=10):

    emissions, global_probs, marginal_probs =  get_averaged_emissions_and_probs(model_path, input_ids, input_mask, batch_size)

    viterbi_input_mask = make_viterbi_mask(input_mask)
    viterbi_paths = get_viterbi_paths(model_path, emissions, viterbi_input_mask)

    # to numpy
    return global_probs.detach().cpu().numpy(), marginal_probs.detach().cpu().numpy(), viterbi_paths.detach().cpu().numpy()


