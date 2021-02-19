import itertools
bases = ['A', 'C', 'G', 'T']

#Neighborhood generation for any size d for mismatches. Not needed usually 
def generate_neighborhood(kmer, seq, i, n, neighborhood):
    print(i, n, seq)
    if i > len(kmer):
        pass
    elif n == 0:
        tmp = seq + kmer[i:]
        #print(tmp)
        neighborhood[tmp] = True
    else:
        for j in range(i, len(kmer)):
            #print(j)
            for base in bases:
                if kmer[j] != base:
                    tmp = seq + kmer[i: j] + base
                    generate_neighborhood(kmer, tmp, j + 1, n - 1, neighborhood)


def generate_neighborhood_d(kmer, d):
    neighborhood = {}
    generate_neighborhood(kmer, '', 0, d, neighborhood)
    print(len(neighborhood))





#Brute force approach for Implanted MOTIF PROBLEM: Find all (k, d)-motifs in a collection of strings.

def combination(k):
    return (''.join(p) for p in itertools.product('ATCG', repeat=k))

def hamming_distance(pattern, seq):
    return sum(c1 != c2 for c1, c2 in zip(pattern, seq))

def window(s, k):
    for i in range(1 + len(s) - k):
        yield s[i:i+k]

def motif_enumeration(k, d, DNA):
    pattern = set()
    for combo in combination(k):
        if all(any(hamming_distance(combo, pat) <= d
                for pat in window(string, k)) for string in DNA):
            pattern.add(combo)
    return pattern

#Use lines below to use
# dna = """ACGT
# ACGT
# ACGT"""

# dna = dna.split("\n")
# print(motif_enumeration(3, 0, dna))






# A brute force algorithm for the Motif Finding Problem (referred to as BruteForceMotifSearch) 
# considers every possible choice of k-mers Motifs from Dna (one k-mer from each string of n nucleotides) 
# and returns the collection Motifs having minimum score. Because there are n - k + 1 choices of k-mers in each 
# of t sequences, there are (n - k + 1)t different ways to form Motifs. For each choice of Motifs, the algorithm 
# calculates Score(Motifs), which requires k · t steps. Thus, assuming that k is smaller than n, the overall running 
# time of the algorithm is O(nt · k · t). We need to come up with a faster algorithm!




#Now we will have median string problem
#To see why we reformulated the Motif Finding Problem as the equivalent Median String Problem, 
# consider the runtime of MedianString and BruteForceMotifSearch. 
# The former algorithm computes d(Pattern, Dna) for each of the 4k k-mers Pattern. 
# Each computation of d(Pattern, Dna) requires a single pass over each string in Dna, 
# which requires approximately k · n · t operations for t strings of length n in Dna. 
# Therefore, MedianString has a running time of O(4^k · n · k · t), 
# which in practice compares favorably with the O(n^t · k · t) running time of BruteForceMotifSearch because the 
# length of a motif (k) typically does not exceed 20 nucleotides, whereas t is measured in the thousands.
#Of course, the ultimate test of a bioinformatics algorithm is how it performs in practice. 
# Unfortunately, since MedianString has to consider 4k k-mers, it becomes too slow for the Subtle Motif Problem, 
# for which k = 15. We will run MedianString with k = 13 in the hope that it will capture a 
# substring of the correct 15-mer motif. 


def generate_kmers(kmer, k, kmers):
    if len(kmer) == k:
        kmers[kmer] = len(kmer) + 1
    else:
        for base in bases:
            generate_kmers(kmer + base, k, kmers)

def generate_kmers_k(k):
    kmer = ''
    kmers = {}
    for base in bases:
        generate_kmers(kmer + base, k, kmers)
    return kmers

def median_string(k, dnas):
    kmers = generate_kmers_k(k)
    print('We have', len(kmers), ' kmer candidates for the median.')
    for kmer in kmers:
        kmers[kmer] = [k + 1] * len(dnas)
    for d, dna in enumerate(dnas): # for each DNA sequence
        for i in range(len(dna) - k + 1): # enumerate kmers
            kmer = dna[i: i + k]
            for pattern in kmers: # calculate hamming distance with each pattern
                kmers[pattern][d] = min(kmers[pattern][d], hamming_distance(pattern, kmer))
    mean_kmer = min(kmers, key = lambda kmer: sum(kmers[kmer]))
    # check if there are multiple answers
    d = sum(kmers[mean_kmer])
    for kmer in kmers:
        if sum(kmers[kmer]) == d:
            print(kmer)    


#Usage of median string problem
# dnas = """AAATTGACGCAT
# GACGACCACGTT
# CGTCAGCGCCTG
# GCTGAGCACCGG
# AGTACGGGACAG"""

# dnas = dnas.split("\n")
# median_string(3, dnas)









#GREEDY MOTIF SEARCH - ORIGINAL


def ProbableKmer(string, matrix):
    probable = 1
    for i in range(len(string)):
        if string[i] == 'A':
            probable *= matrix[0][i]
        if string[i] == 'C':
            probable *= matrix[1][i]
        if string[i] == 'G':
            probable *= matrix[2][i]
        if string[i] == 'T':
            probable *= matrix[3][i]
    return probable

# Profile-most probable k-mer in the i-th string in Dna
def FindProfileMostProbableKmer(string, k, matrix):
    seq = {}
    for i in range(len(string) - k + 1):
        seq[string[i:i + k]] = ProbableKmer(string[i:i + k], matrix)
    max_key = sorted(seq.items(), key=lambda x:x[1], reverse=True)[0][0]
    return max_key

# Score(Motifs)
def Score(Motifs):
    score = 0
    for i in range(len(Motifs[0])):
        j = [motif[i] for motif in Motifs]
        score += (len(j) - max(j.count("A"), j.count("C"), j.count("T"), j.count("G")))
    return score

def GreedyMotifSearch(Dna, k, t):
    # BestMotifs ← motif matrix formed by first k-mers in each string from Dna
    BestMotifs = [dna[:k] for dna in Dna]
    # for each k-mer Motif in the first string from Dna
    for k_mer in [Dna[0][i:i+k] for i in range(len(Dna[0])-k+1)]:
        # Motif1 ← Motif
        Motifs = [k_mer]
        # for i = 2 to t
        for i in range(1, t):
            # form Profile from motifs Motif1, …, Motifi - 1
            motifs = Motifs[:i]
            # Motifi ← Profile-most probable k-mer in the i-th string in Dna
            matrix = []
            for nar in ["A", "C", "G", "T"]:
                mat = []
                for j in range(k):
                    mm = [m[j] for m in motifs]
                    mat.append(mm.count(nar)/len(motifs))
                matrix.append(mat)
            # Motifs ← (Motif1, …, Motift)    
            Motifs.append(FindProfileMostProbableKmer(Dna[i], k, matrix))
        # print(Motifs)
        # if Score(Motifs) < Score(BestMotifs), BestMotifs ← Motifs
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs




#usage for greedy motif 
# dnas = """AAATTGACGCAT
# GACGACCACGTT
# CGTCAGCGCCTG
# GCTGAGCACCGG
# AGTACGGGACAG"""

# dnas = dnas.split("\n")

#print(GreedyMotifSearch(dnas, 3, 5))




#usage for Find a Profile-most Probable k-mer in a String
# matrix = """0.2 0.2 0.3 0.2 0.3
# 0.4 0.3 0.1 0.5 0.1
# 0.3 0.3 0.5 0.2 0.4
# 0.1 0.2 0.1 0.1 0.2"""

# matrix = matrix.split("\n")
# matrix_f = [[float(l) for l in x.strip().split()] for x in matrix]

# print(FindProfileMostProbableKmer("ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT", 5 , matrix_f))








#GREEDY MOTIF with psuedocounts


def profile_most_probable_kmer(dna, k, profile):
    '''Returns the profile most probable k-mer for the given input data.'''
    # A dictionary relating nucleotides to their position within the profile.
    nuc_loc = {nucleotide:index for index,nucleotide in enumerate('ACGT')}

    # Initialize the maximum probabily.
    max_probability = -1

    # Compute the probability of the each k-mer, store it if it's currently a maximum.
    for i in range(len(dna)-k+1):
        # Get the current probability.
        current_probability = 1
        for j, nucleotide in enumerate(dna[i:i+k]):
            current_probability *= profile[j][nuc_loc[nucleotide]]

        # Check for a maximum.
        if current_probability > max_probability:
            max_probability = current_probability
            most_probable = dna[i:i+k]

    return most_probable


def score(motifs):
    '''Returns the score of the given list of motifs.'''
    columns = [''.join(seq) for seq in zip(*motifs)]
    max_count = sum([max([c.count(nucleotide) for nucleotide in 'ACGT']) for c in columns])
    return len(motifs[0])*len(motifs) - max_count


def profile_with_pseudocounts(motifs):
    '''Returns the profile of the dna list motifs.'''
    columns = [''.join(seq) for seq in zip(*motifs)]
    return [[float(col.count(nuc)+1) / float(len(col)+4) for nuc in 'ACGT'] for col in columns]


def greedy_motif_search_pseudocounts(dna_list, k, t):
    '''Runs the Greedy Motif Search with Pseudocounts algorithm and retuns the best motif.'''
    # Initialize the best score as a score higher than the highest possible score.
    best_score = t*k

    # Run the greedy motif search.
    for i in range(len(dna_list[0]) - k + 1):
        # Initialize the motifs as each k-mer from the first dna sequence.
        motifs = [dna_list[0][i:i + k]]

        # Find the most probable k-mer in the next string, using pseudocounts.
        for j in range(1, t):
            current_profile = profile_with_pseudocounts(motifs)
            motifs.append(profile_most_probable_kmer(dna_list[j], k, current_profile))

        # Check to see if we have a new best scoring list of motifs.
        current_score = score(motifs)
        if current_score < best_score:
            best_score = current_score
            best_motifs = motifs

    # for x in best_motifs:
    #     print(x)

    return best_motifs        


#usage for greedy motif  WITH PSUEDOCOUNTS
# dnas = """GGCGTTCAGGCA
# AAGAATCAGTCA
# CAAGGAGTTCGC
# CACGTCAATCAC
# CAATAATATTCG"""
# dnas = dnas.split("\n")

# print(greedy_motif_search_pseudocounts(dnas, 3, 5))

