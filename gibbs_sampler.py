# Implement the Gibbs sampler for motif finding in DNA sequences

import sys
sys.path.append('Arthur code')
from arthur_gibbs_sampler import gibbs_sampler
from arthur_greedy_motif_search_pseudocounts import score_motifs_pc

def read_fasta(filename: str) -> str:
    """Read a FASTA file and return the concatenated DNA sequence as a string.
    
    Args:
    filename (str): The path to the FASTA file.

    Returns:
        sequence (str): The concatenated DNA sequence.
    """
    
    sequence = ""
    with open(filename, "r") as file:
        for line in file:
            if line.startswith(">"):
                continue  # Skip header lines
            else:
                sequence += line.strip()  # Remove whitespace and concatenate
            
    return sequence


def hamming_distance(seq1: str, seq2: str) -> int:
    """Calculate the Hamming distance, number of positions at which the corresponding symbols are different, between two sequences.
    
    Args:
        seq1 (str): The first sequence.
        seq2 (str): The second sequence.

    Returns:
        distance (int): The Hamming distance between the two sequences.
    """
    
    return sum(el1 != el2 for el1, el2 in zip(seq1, seq2))


def gibbs_sampler_per_operon(operon_groups: dict, k: int = 10, t: int = 6, num_runs: int = 100) -> dict:
    """
    Run the Gibbs sampler on a given FASTA file and return the identified motifs.

    Args:
        operon_groups (dict): A dictionary mapping operon names to lists of FASTA files.
        k (int): The length of the motif to find.
        t (int): The number of sequences to consider.
        num_runs (int): The number of iterations for the Gibbs sampler.

    Returns:
        all_results (dict): A dictionary mapping operon names to lists of identified motifs.
    """
    
    all_results = {}
    
    # Go through each operon and its associated FASTA files
    for operon, filenames in operon_groups.items():
        # Help from Claude Code Python-ifying this for-loop
        sequences = [read_fasta(filename) for filename in filenames]
        best_motif = None
        best_score = float('inf')
        
        # Run the Gibbs sampler multiple times to find the best motif
        for run in range(num_runs):
            motifs = gibbs_sampler(sequences, k, t, num_runs)
            score = score_motifs_pc(motifs)
            if score < best_score:
                best_score = score
                best_motif = motifs[:]
    
        all_results[operon] = {
                "motif": best_motif,
                "score": best_score,
                "sequences": filenames
            }
    
    return all_results



