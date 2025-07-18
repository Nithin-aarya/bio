# 1. Count the ATGC content of a given DNA sequence
def count_nucleotides(dna_sequence):
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for nucleotide in dna_sequence:
        counts[nucleotide] += 1
    return counts

# 2. Return the complement of a DNA strand
def complement_dna(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement[base] for base in dna_sequence])

# 3. Convert a DNA sequence into an RNA sequence
def dna_to_rna(dna_sequence):
    return dna_sequence.replace('T', 'U')

# 4. Translate RNA into protein sequence
def translate_rna(rna_sequence):
    genetic_code = {
        'AUG': 'Methionine', 'UUU': 'Phenylalanine', 'UUC': 'Phenylalanine',
        # Add remaining codons as needed
    }
    return [genetic_code[rna_sequence[i:i+3]] for i in range(0, len(rna_sequence), 3)]

# 5. Find all occurrences of a motif in a DNA sequence
def find_motif(dna_sequence, motif):
    positions = []
    start = 0
    while True:
        start = dna_sequence.find(motif, start) + 1
        if start > 0:
            positions.append(start)
        else:
            break
    return positions

# 6. Calculate the total mass of a protein sequence
def calculate_protein_mass(protein_sequence):
    amino_acid_mass = {
        'A': 71.03711, 'C': 103.00919,
        # Add other amino acids
    }
    return sum([amino_acid_mass[aa] for aa in protein_sequence])

# 7. Count the number of occurrences of each nucleotide
# Same as #1

# 8. Find most frequent k-mers in a DNA sequence
def find_most_frequent_kmers(dna_sequence, k):
    kmers_count = {}
    for i in range(len(dna_sequence) - k + 1):
        kmer = dna_sequence[i:i+k]
        if kmer in kmers_count:
            kmers_count[kmer] += 1
        else:
            kmers_count[kmer] = 1
    max_count = max(kmers_count.values())
    most_frequent = [kmer for kmer, count in kmers_count.items() if count == max_count]
    return most_frequent, max_count

# 9. Find reverse complement of a DNA Sequence
def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([complement[n] for n in reversed(dna_sequence)])

# 10. Design primers
def calculate_tm(seq):
    a_t = seq.count('A') + seq.count('T')
    g_c = seq.count('G') + seq.count('C')
    return 2 * a_t + 4 * g_c

def gc_content(seq):
    g_c = seq.count('G') + seq.count('C')
    return (g_c / len(seq)) * 100

def design_primer(dna_sequence, primer_length=20):
    primers = []
    for i in range(len(dna_sequence) - primer_length + 1):
        primer = dna_sequence[i:i+primer_length]
        if 40 <= gc_content(primer) <= 60:
            primers.append(primer)
    return primers

# 11. Find pattern (e.g., origin of replication)
def find_pattern_in_dna(sequence, pattern):
    indices = []
    for i in range(len(sequence) - len(pattern) + 1):
        if sequence[i:i+len(pattern)] == pattern:
            indices.append(i)
    return indices

# 12. Evaluate Open Reading Frames (ORFs)
def find_orfs(dna_seq):
    def reverse_complement(dna_seq):
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return ''.join([complement[base] for base in reversed(dna_seq)])

    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    for strand, seq in [('+', dna_seq), ('-', reverse_complement(dna_seq))]:
        for frame in range(3):
            trans_start_pos = None
            for pos in range(frame, len(seq) - 2, 3):
                codon = seq[pos:pos + 3]
                if trans_start_pos is None and codon == start_codon:
                    trans_start_pos = pos
                elif trans_start_pos is not None and codon in stop_codons:
                    orfs.append((strand, frame, trans_start_pos, pos + 3))
                    trans_start_pos = None
    return orfs

# 13. Count RFLP markers
def count_rflp_markers(dna_sequence, restriction_enzymes):
    rflp_marker_counts = {}
    for enzyme, cut_site in restriction_enzymes.items():
        rflp_marker_counts[enzyme] = dna_sequence.count(cut_site)
    return rflp_marker_counts

# 14. Evaluate mutation frequency
def evaluate_mutations(sequence1, sequence2):
    if len(sequence1) != len(sequence2):
        raise ValueError("Sequences must be of equal length.")
    mutation_count = sum(1 for b1, b2 in zip(sequence1, sequence2) if b1 != b2)
    mutation_frequency = (mutation_count / len(sequence1)) * 100
    return mutation_count, mutation_frequency

# 15. Calculate mass of the peptide
def calculate_peptide_mass(sequence):
    mass_table = {
        'A': 71.0788, 'R': 156.1875, 'N': 114.1038, 'D': 115.0886,
        'C': 103.1388, 'E': 129.1155, 'Q': 128.1307, 'G': 57.0519,
        'H': 137.1411, 'I': 113.1594, 'L': 113.1594, 'K': 128.1741,
        'M': 131.1926, 'F': 147.1766, 'P': 97.1167, 'S': 87.0782,
        'T': 101.1051, 'W': 186.2132, 'Y': 163.1760, 'V': 99.1326
    }
    water_mass = 18.015
    total_mass = sum(mass_table[aa] for aa in sequence)
    total_mass -= (len(sequence) - 1) * water_mass
    return total_mass

# 16. Dot plot
import matplotlib.pyplot as plt
def create_dot_plot(seq1, seq2):
    matrix = [[1 if seq1[i] == seq2[j] else 0 for j in range(len(seq2))] for i in range(len(seq1))]
    plt.imshow(matrix, cmap='Greys', interpolation='none')
    plt.title('Dot Plot')
    plt.xlabel('Sequence 2')
    plt.ylabel('Sequence 1')
    plt.show()

# 17. Hamming distance
def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))

# 18. Levenshtein distance
def levenshtein_distance(seq1, seq2):
    if len(seq1) < len(seq2):
        return levenshtein_distance(seq2, seq1)
    if len(seq2) == 0:
        return len(seq1)
    previous_row = range(len(seq2) + 1)
    for i, c1 in enumerate(seq1):
        current_row = [i + 1]
        for j, c2 in enumerate(seq2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]

# 19. Global alignment (Needleman-Wunsch)
def needleman_wunsch(seq1, seq2, match_score=1, mismatch_score=-1, gap_penalty=-2):
    m, n = len(seq1), len(seq2)
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    traceback_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(n + 1):
        score_matrix[0][j] = j * gap_penalty
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_score)
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            score_matrix[i][j] = max(match, delete, insert)
            if score_matrix[i][j] == match:
                traceback_matrix[i][j] = "diag"
            elif score_matrix[i][j] == delete:
                traceback_matrix[i][j] = "up"
            else:
                traceback_matrix[i][j] = "left"
    align1, align2 = "", ""
    i, j = m, n
    while i > 0 or j > 0:
        if traceback_matrix[i][j] == "diag":
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == "up":
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1
    return align1, align2

# 20. Local alignment (Smith-Waterman)
def smith_waterman(seq1, seq2, match_score=2, gap_cost=1):
    m, n = len(seq1), len(seq2)
    score = [[0 for _ in range(n+1)] for _ in range(m+1)]
    traceback = [[0 for _ in range(n+1)] for _ in range(m+1)]
    max_score = 0
    max_pos = None
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = score[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else -gap_cost)
            delete = score[i-1][j] - gap_cost
            insert = score[i][j-1] - gap_cost
            score[i][j] = max(0, match, delete, insert)
            if score[i][j] == match:
                traceback[i][j] = '↖'
            elif score[i][j] == delete:
                traceback[i][j] = '↑'
            elif score[i][j] == insert:
                traceback[i][j] = '←'
            else:
                traceback[i][j] = None
            if score[i][j] >= max_score:
                max_score = score[i][j]
                max_pos = (i, j)
    align1, align2 = '', ''
    i, j = max_pos
    while traceback[i][j] is not None:
        if traceback[i][j] == '↖':
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1
        elif traceback[i][j] == '↑':
            align1 = seq1[i-1] + align1
            align2 = '-' + align2
            i -= 1
        elif traceback[i][j] == '←':
            align1 = '-' + align1
            align2 = seq2[j-1] + align2
            j -= 1
    return align1, align2, max_score
