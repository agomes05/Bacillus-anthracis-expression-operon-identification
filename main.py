# 1. Write read_fasta() to return a single sequence string
# 2. Create group_files_by_operon() to organize your 18 files into 3 groups
# 3. Write run_gibbs_by_operon() to run Gibbs on each group (6 sequences each)
# 4. Write hamming_distance() and clustering code
# 5. Compare the 3 resulting motifs to answer: "Are gerA, gerB, gerK regulated the same way?"
# 6. Validate: Check that each motif appears in both species

from gibbs_sampler import gibbs_sampler

# --- CONSTANTS --- #

# Divided by operons because we want to see the consistency across species for each operon
FILES_TO_PROCESS = {
        "gerA": [
            "BA-GERA250BP.fasta",
            "BA-GERA500BP.fasta",
            "BA-GERA1000BP.fasta",
            "BC-GERA250BP.fasta",
            "BC-GERA500BP.fasta",
            "BC-GERA1000BP.fasta",
        ],
        "gerB": [
            "BA-GERB250BP.fasta",
            "BA-GERB500BP.fasta",
            "BA-GERB1000BP.fasta",
            "BC-GERB250BP.fasta",
            "BC-GERB500BP.fasta",
            "BC-GERB1000BP.fasta",
        ],
        "gerK": [
            "BA-GERK250BP.fasta",
            "BA-GERK500BP.fasta",
            "BA-GERK1000BP.fasta",
            "BC-GERK250BP.fasta",
            "BC-GERK500BP.fasta",
            "BC-GERK1000BP.fasta",
        ]
    }

K = 10  # Length of the motif to find

if __name__ == "__main__":
    for filename, species, operon, sequence_size in FILES_TO_PROCESS:
        print(f"Processing {filename} for {species} {operon} with sequence size {sequence_size}...")
        results = gibbs_sampler(filename, species, operon, sequence_size)
        print(f"Results for {filename}: {results}\n")