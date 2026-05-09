# 1. Write read_fasta() to return a single sequence string
# 2. Create group_files_by_operon() to organize your 18 files into 3 groups
# 3. Write run_gibbs_by_operon() to run Gibbs on each group (6 sequences each)
# 4. Write hamming_distance() and clustering code
# 5. Compare the 3 resulting motifs to answer: "Are gerA, gerB, gerK regulated the same way?"
# 6. Validate: Check that each motif appears in both species

import os
from gibbs_sampler import gibbs_sampler_per_operon

base_dir = os.path.dirname(__file__)
anthracis_dir = os.path.join(base_dir, "B. anthracis sequences")
cereus_dir = os.path.join(base_dir, "B. cereus sequences")

# --- CONSTANTS --- #

# Divided by operons because we want to see the consistency across species for each operon
FILES_TO_PROCESS = {
    "gerA": [
        os.path.join(anthracis_dir, "BA-GERA250BP.fasta"),
        os.path.join(anthracis_dir, "BA-GERA500BP.fasta"),
        os.path.join(anthracis_dir, "BA-GERA1000BP.fasta"),
        os.path.join(cereus_dir, "BC-GERA250BP.fasta"),
        os.path.join(cereus_dir, "BC-GERA500BP.fasta"),
        os.path.join(cereus_dir, "BC-GERA1000BP.fasta"),
    ],
    "gerB": [
        os.path.join(anthracis_dir, "BA-GERB250BP.fasta"),
        os.path.join(anthracis_dir, "BA-GERB500BP.fasta"),
        os.path.join(anthracis_dir, "BA-GERB1000BP.fasta"),
        os.path.join(cereus_dir, "BC-GERB250BP.fasta"),
        os.path.join(cereus_dir, "BC-GERB500BP.fasta"),
        os.path.join(cereus_dir, "BC-GERB1000BP.fasta"),
    ],
    "gerK": [
        os.path.join(anthracis_dir, "BA-GERK250BP.fasta"),
        os.path.join(anthracis_dir, "BA-GERK500BP.fasta"),
        os.path.join(anthracis_dir, "BA-GERK1000BP.fasta"),
        os.path.join(cereus_dir, "BC-GERK250BP.fasta"),
        os.path.join(cereus_dir, "BC-GERK500BP.fasta"),
        os.path.join(cereus_dir, "BC-GERK1000BP.fasta"),
    ],
}

K = 10  # Length of the motif to find

if __name__ == "__main__":
    for operon, filenames in FILES_TO_PROCESS.items():
        for filename in filenames:
            print(f"Processing {filename} for {operon}...")
            results = gibbs_sampler_per_operon({operon: [filename]}, K, t=6, num_runs=100)
            print(f"Results for {filename}: {results}\n")