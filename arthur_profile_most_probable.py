import numpy as np

DNA = "ACGT"

def kmer_probability(kmer, profile):
    """Return the probability of a k-mer given a profile matrix"""
    # Start with largest possibility
    prob = 1.0
    # Multiply probabilities of each nucleotide in k-mer based on profile matrix
    for i, nucleotide in enumerate(kmer):
        prob *= profile[DNA.index(nucleotide), i]
        
    return prob

def profile_most_probable_kmer(sequence, k, profile):
    """Return the Profile-most probable k-mer in sequence"""
    # Initialize probability to keep track of the most probable probability
    best_prob = -1.0
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        prob = kmer_probability(kmer, profile)
        if prob > best_prob:
            best_prob = prob
            most_probable_kmer = kmer
            
    return most_probable_kmer
    
def read2C(filename):
    """Read 2C input data"""
    with open(filename, "r") as file:
        sequence = file.readline().strip()
        k = int(file.readline().strip())
        profile = np.loadtxt(file)  # Default type is float
        
        return sequence, k, profile

# Usage
if __name__ == "__main__":
    sequence, k, profile = read2C("/Users/arthurgomes/Documents/CSCI/CSCI321/Rosalind 1/Profile-most Probable k-mer/rosalind_ba2c.txt")
    result = profile_most_probable_kmer(sequence, k, profile)
    print(result)