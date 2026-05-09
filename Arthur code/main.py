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

if __name__ == "__main__":
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