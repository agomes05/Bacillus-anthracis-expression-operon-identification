import sys
sys.path.append('/Users/arthurgomes/Documents/CSCI/CSCI321/Rosalind 1/GreedyMotifSearch with Pseudocounts')
from arthur_greedy_motif_search_pseudocounts import build_profile_pc, score_motifs_pc
sys.path.append('/Users/arthurgomes/Documents/CSCI/CSCI321/Rosalind 1/Profile-most Probable k-mer')
from arthur_profile_most_probable import profile_most_probable_kmer
import random

def gibbs_sampler(dna, k, t, N):
    """Gibbs Sampler algorithm"""
    # Randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
    # Initialize BestMotifs ← Motifs
    motifs = []
    for seq in dna:
        start = random.randint(0, len(seq) - k)
        motifs.append(seq[start:start + k])
    best_motifs = motifs[:]
    
    for j in range(N):
        i = random.randint(0, t-1)
        # Profile ← profile matrix constructed from all strings in Motifs except for Motifi
        profile = build_profile_pc(motifs[:i] + motifs[i+1:])
        motifi = profile_most_probable_kmer(dna[i], k, profile)
        # Motifi is the most probable k-mer in the i-th string according to Profile
        motifs[i] = motifi
        if score_motifs_pc(motifs) < score_motifs_pc(best_motifs):
            best_motifs = motifs[:]
    return best_motifs

# Debugged with the help of Claude Code
def read2D(filename):
    with open(filename, "r") as file:
        data = file.read().split()
    k = int(data[0])
    t = int(data[1])
    N = int(data[2])
    dna = data[3:]
    return k, t, N, dna

if __name__ == "__main__":
    k, t, N, dna = read2D("/Users/arthurgomes/Documents/CSCI/CSCI321/Rosalind 2/Gibbs Sampler/GibbsSampler/inputs/input_1.txt")
    
    best_motifs = None
    best_score = float('inf')
    
    # A collection BestMotifs resulting from running GibbsSampler(Dna, k, t, N) 20 times.
    for i in range(20):
        motifs = gibbs_sampler(dna, k, t, N)
        score = score_motifs_pc(motifs)
        if score < best_score:
            best_score = score
            best_motifs = motifs[:]
    
    print("\n".join(best_motifs))