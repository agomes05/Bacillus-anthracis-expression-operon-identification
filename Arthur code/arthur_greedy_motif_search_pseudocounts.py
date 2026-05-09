import numpy as np

# Make sure this path still points to your profile_most_probable function!
from arthur_profile_most_probable import profile_most_probable_kmer

DNA = "ACGT"

def build_profile_pc(motifs):
    """Build a profile matrix from a list of motifs (with pseudocounts)"""
    profile = np.ones((4, len(motifs[0])), dtype=float)
    for motif in motifs:
        for i, nucleotide in enumerate(motif):
            profile[DNA.index(nucleotide), i] += 1
    # Add 4 to account for pseudocounts (one per nucleotide)
    profile /= (len(motifs) + 4)
    return profile
    
def score_motifs_pc(motifs):
    """Return the score of a motif matrix using pseudocount profiles"""
    # MUST call the pseudocount version of build_profile
    profile = build_profile_pc(motifs)
    consensus = "".join(DNA[np.argmax(profile[:, i])] for i in range(profile.shape[1]))
    score = sum(sum(1 for motif in motifs if motif[i] != consensus[i]) for i in range(len(consensus)))
    return score

def greedy_motif_search_pc(dna, k, t):
    """GreedyMotifSearch using Pseudocounts"""
    best_motifs = [seq[:k] for seq in dna]
    for i in range(len(dna[0]) - k + 1):
        motif1 = dna[0][i:i+k]
        motifs = [motif1]
        
        for j in range(1, len(dna)):
            # MUST call the pseudocount version here too
            profile = build_profile_pc(motifs)
            motif_j = profile_most_probable_kmer(dna[j], k, profile)
            motifs.append(motif_j)
            
        # MUST use the pseudocount scoring function
        if score_motifs_pc(motifs) < score_motifs_pc(best_motifs):
            best_motifs = motifs
            
    return best_motifs
    
def read2E(filename):
    """Read rosalind.info input data safely, handling all whitespace"""
    with open(filename, "r") as file:
        data = file.read().split()
        k = int(data[0])
        t = int(data[1])
        dna = data[2:] 
    return k, t, dna

if __name__ == "__main__":
    k, t, dna = read2E("/Users/arthurgomes/Documents/CSCI/CSCI321/Rosalind 1/GreedyMotifSearch with Pseudocounts/rosalind_ba2e.txt")
    
    best_motifs = greedy_motif_search_pc(dna, k, t)
    print("\n".join(best_motifs))