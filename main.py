# 1. Run Gibbs Sampler on all 18 sequence files (or a subset) → get ~18 motifs
# 2. Compute distances between all pairs of motifs (how different are they?)
# 3. Run hierarchical clustering on those distances
# 4. Identify the largest cluster → this is your conserved motif
# 5. Compare cluster membership across species:

# If B. anthracis and B. cereus motifs are in the SAME cluster → conserved!
# If they're in DIFFERENT clusters → diverged!

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