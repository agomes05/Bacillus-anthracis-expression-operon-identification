from tqdm import tqdm
from gibbs_sampler import gibbs_sampler_per_operon
from hierarchichal_clustering import hierarchical_clustering
import pandas as pd
from gibbs_sampler import distance_matrix, gibbs_sampler_per_operon

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
    all_results = gibbs_sampler_per_operon(FILES_TO_PROCESS, k=10, t=6, num_runs=2000, allotted_memory = 150)
    motif_rows = []
    for operon, results in all_results.items():
        print(f"\nResults for {operon}:")
        for filename, motif in zip(results["sequences"], results["motif"]):
            print(filename, motif)
            motif_rows.append({
                "operon": operon, "filename": filename, "motif": motif, "score": results["score"]})
        print("Score:", results["score"])
        matrix = matrix_for_clustering(results["motif"])
        CSV_formatrix = pd.DataFrame(matrix, index=results["sequences"], columns=results["sequences"])
        CSV_formatrix.to_csv(f"{operon}_MATRIX_OUTPUT_FOR_DATA_AC_AG.csv")
        print("distance matrix:")
        for row in matrix:
            print(row)
        clusters = hierarchical_clustering(len(results["motif"]), matrix)
        print("clusters:")
        for cluster in clusters:
            print(cluster)
    motif_dataframe = pd.DataFrame(motif_rows)
    motif_dataframe.to_csv("FinalMotifs.csv", index = False)
    representative_motifs = {
        "gerA_BA": all_results["gerA"]["motif"][0],
        "gerA_BC": all_results["gerA"]["motif"][3],
        "gerB_BA": all_results["gerB"]["motif"][0],
        "gerB_BC": all_results["gerB"]["motif"][3],
        "gerK_BA": all_results["gerK"]["motif"][0],
        "gerK_BC": all_results["gerK"]["motif"][3],}
    labels = list(representative_motifs.keys())
    motif_values = list(representative_motifs.values())
    dist_matrix = matrix_for_clustering(motif_values)
    dist_matrix_for_df = pd.DataFrame(dist_matrix, index=labels, columns=labels)
    dist_matrix_for_df.to_csv("Distance_matrix_representation.csv")
    print("\nRepresentative motif distance matrix:")
    print(dist_matrix_for_df)
    clusters = hierarchical_clustering(len(motif_values), dist_matrix)
    print("\nRepresentative motif clustering:")
    for cluster in clusters:
        print([labels[i - 1] for i in cluster])
