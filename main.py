# 1. Run Gibbs Sampler on all 18 sequence files (or a subset) → get ~18 motifs
# 2. Compute distances between all pairs of motifs (how different are they?)
# 3. Run hierarchical clustering on those distances
# 4. Identify the largest cluster → this is your conserved motif
# 5. Compare cluster membership across species:

# If B. anthracis and B. cereus motifs are in the SAME cluster → conserved!
# If they're in DIFFERENT clusters → diverged!