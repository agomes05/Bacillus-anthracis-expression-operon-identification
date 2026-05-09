<<<<<<< HEAD
# 1. Write read_fasta() to return a single sequence string
# 2. Create group_files_by_operon() to organize your 18 files into 3 groups
# 3. Write run_gibbs_by_operon() to run Gibbs on each group (6 sequences each)
# 4. Write hamming_distance() and clustering code
# 5. Compare the 3 resulting motifs to answer: "Are gerA, gerB, gerK regulated the same way?"
# 6. Validate: Check that each motif appears in both species
import sys
sys.path.append("/Users/alexchrostowski/Desktop/Bioinformatics Algorithms ")
 
from gibbs_sampler import gibbs_sampler_per_operon
from HierarchalClustering import Hierarchal_clustering
=======
import pandas as pd
from gibbs_sampler import distance_matrix, gibbs_sampler_per_operon
>>>>>>> d53cd785ec44d518503272de031b2c1665fe0f1b

# --- CONSTANTS --- #

# Divided by operons because we want to see the consistency across species for each operon
FILES_TO_PROCESS = {
    "gerA": [
        "B. anthracis sequences/BA-GERA250BP.fasta",
        "B. anthracis sequences/BA-GERA500BP.fasta",
        "B. anthracis sequences/BA-GERA1000BP.fasta",
        "B. cereus sequences/BC-GERA250BP.fasta",
        "B. cereus sequences/BC-GERA500BP.fasta",
        "B. cereus sequences/BC-GERA1000BP.fasta",
    ],
    "gerB": [
        "B. anthracis sequences/BA-GERB250BP.fasta",
        "B. anthracis sequences/BA-GERB500BP.fasta",
        "B. anthracis sequences/BA-GERB1000BP.fasta",
        "B. cereus sequences/BC-GERB250BP.fasta",
        "B. cereus sequences/BC-GERB500BP.fasta",
        "B. cereus sequences/BC-GERB1000BP.fasta",
    ],
    "gerK": [
        "B. anthracis sequences/BA-GERK250BP.fasta",
        "B. anthracis sequences/BA-GERK500BP.fasta",
        "B. anthracis sequences/BA-GERK1000BP.fasta",
        "B. cereus sequences/BC-GERK250BP.fasta",
        "B. cereus sequences/BC-GERK500BP.fasta",
        "B. cereus sequences/BC-GERK1000BP.fasta",
    ]
}

K = 10  # Length of the motif to find

def mod_hamming_distance(motif1, motif2):
    return sum(motif1[i] != motif2[i] for i in range(len(motif1)))
def matrix_for_clustering(motifs):
    matrix = []
    for motif1 in motifs: 
        row = []
        for motif2 in motifs:
            row.append(float(mod_hamming_distance(motif1,motif2)))
        matrix.append(row)
    return matrix 

if __name__ == "__main__":
<<<<<<< HEAD
    all_results = gibbs_sampler_per_operon(FILES_TO_PROCESS, k=10, t=6, num_runs=200)
    for operon, results in all_results.items():
        print(f"\nResults for {operon}:")
        for filename, motif in zip(results["sequences"], results["motif"]):
            print(filename, motif)
        print("Score:", results["score"])
        matrix = matrix_for_clustering(results["motif"])
        print("distance matrix:")
        for row in matrix: 
            print(row)
        clusters = Hierarchal_clustering(len(results["motif"]), matrix)
        print("clusters:")
        for cluster in clusters:
            print(cluster)
=======
    all_results = gibbs_sampler_per_operon(FILES_TO_PROCESS, k=10, t=6, num_runs=20)
    motifs = {
        "gerA_BA": all_results["gerA"]["motif"][0], # Indices 0-2 are the motifs from B. anthracis
        "gerA_BC": all_results["gerA"]["motif"][3], # Indices 3-5 are the motifs from B. cereus
        "gerB_BA": all_results["gerB"]["motif"][0],
        "gerB_BC": all_results["gerB"]["motif"][3],
        "gerK_BA": all_results["gerK"]["motif"][0],
        "gerK_BC": all_results["gerK"]["motif"][3]
    }
    print("Motifs identified for each operon and species:")
    for key, motif in motifs.items():
        print(f"{key}: {motif}")
    
    dist_matrix = distance_matrix(motifs)
    print("\nDistance Matrix:")
    dist_matrix_df = pd.DataFrame(dist_matrix, index=motifs.keys(), columns=motifs.keys())
    print(dist_matrix_df)
>>>>>>> d53cd785ec44d518503272de031b2c1665fe0f1b
