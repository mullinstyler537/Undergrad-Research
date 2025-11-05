"""
This script calculates omega (dN/dS) from a two-sequence alignment
The first sequence in the alignment is used as the reference sequence
Exon regions must be provided as intervals in a csv file 
  - the first column is the first position of an exon and the second column shows the final position of the exon

--------------

Example usage:

python3 calculate-omega.py selection-alignment.fasta intervals.csv

"""
import csv
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from Bio.Data import CodonTable

# Load codon table for standard translation and mutation comparisons
standard_table = CodonTable.unambiguous_dna_by_id[1]

def read_intervals(interval_file):
    """
    Reads exon intervals from a CSV file.

    Args:
    interval_file (str): The path to the CSV file containing exon intervals. 
                         Each row in the file should have two columns:
                         the start and end positions of the exon.

    Returns:
    list of tuple: A list of intervals where each interval is a tuple (start, end).
                   Positions are converted to 0-based indexing.
    """
    intervals = []
    with open(interval_file, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            # Convert to 0-based indexing (subtract 1 from start position)
            start, end = int(row[0]) - 1, int(row[1])
            intervals.append((start, end))
    return intervals

def splice_sequence(seq, intervals):
    """
    Extracts coding regions from a sequence based on provided exon intervals.

    Args:
    seq (str): The original DNA sequence.
    intervals (list of tuple): List of exon intervals (start, end).

    Returns:
    str: A string representing the spliced sequence (only coding regions).
    """
    # Join the sequence fragments from the specified intervals
    return ''.join([seq[start:end] for start, end in intervals])

def calculate_np_sp(ref_coding_seq):
    """
    Calculates the number of possible nonsynonymous (NP) and synonymous (SP) sites
    in the reference coding sequence.

    Args:
    ref_coding_seq (str): The reference coding DNA sequence.

    Returns:
    tuple: A tuple (NP, SP) where NP is the number of possible nonsynonymous sites
           and SP is the number of possible synonymous sites.
    """
    NP = 0  # Non-synonymous possible sites
    SP = 0  # Synonymous possible sites

    # Loop through the coding sequence codon by codon (3 nucleotides at a time)
    for i in range(0, len(ref_coding_seq) - 2, 3):
        codon = ref_coding_seq[i:i+3]
        aa_ref = Seq(codon).translate(table=standard_table)  # Translate the codon to amino acid
        
        # For each nucleotide in the codon, test possible substitutions
        for pos in range(3):
            original_base = codon[pos]
            for base in 'ACGT':  # Substitute each nucleotide
                if base != original_base:
                    # Create the mutated codon
                    mutated_codon = codon[:pos] + base + codon[pos+1:]
                    aa_mut = Seq(mutated_codon).translate(table=standard_table)
                    # Check if the mutation causes a synonymous or nonsynonymous change
                    if aa_mut == aa_ref:
                        NP += 1/3  # Synonymous site
                    else:
                        SP += 1/3  # Nonsynonymous site

                        
    return NP, SP

def calculate_n_s(ref_coding_seq, target_coding_seq):
    """
    Calculates the number of nonsynonymous (N) and synonymous (S) substitutions
    between two coding sequences.

    Args:
    ref_coding_seq (str): The reference coding DNA sequence.
    target_coding_seq (str): The target coding DNA sequence to compare.

    Returns:
    tuple: A tuple (N, S) where N is the number of nonsynonymous substitutions
           and S is the number of synonymous substitutions.
    """
    N = 0  # Nonsynonymous substitutions
    S = 0  # Synonymous substitutions
    
    # Loop through the coding sequences codon by codon
    for i in range(0, len(ref_coding_seq) - 2, 3):
        ref_codon = ref_coding_seq[i:i+3]
        target_codon = target_coding_seq[i:i+3]
        
        # If the codons are different, check if the mutation is nonsynonymous or synonymous
        if ref_codon != target_codon:
            aa_ref = Seq(ref_codon).translate(table=standard_table)
            aa_target = Seq(target_codon).translate(table=standard_table)
            
            if aa_ref == aa_target:
                N += 1  # Synonymous substitution
            else:
                S += 1  # Nonsynonymous substitution

    return N, S

def calculate_omega(fasta_file, interval_file):
    """
    Main function to calculate omega (ω) based on two aligned sequences and exon intervals.

    Args:
    fasta_file (str): The path to the FASTA file containing two aligned sequences.
                      The first sequence is treated as the reference.
    interval_file (str): The path to the CSV file containing exon intervals.
    
    Returns:
    float: The calculated omega (ω), which is the ratio of nonsynonymous to synonymous
           substitution rates (dN/dS).
    """
    # Load sequences from the FASTA file
    with open(fasta_file, 'r') as f:
        sequences = list(SeqIO.parse(f, 'fasta'))
        ref_seq = str(sequences[0].seq)  # First sequence is the reference
        target_seq = str(sequences[1].seq)  # Second sequence is the target
    
    # Load exon intervals from the CSV file
    intervals = read_intervals(interval_file)
    
    # Splice out untranslated and intron regions from both sequences
    ref_coding_seq = splice_sequence(ref_seq, intervals)
    target_coding_seq = splice_sequence(target_seq, intervals)
    
    # Calculate NP (possible nonsynonymous sites) and SP (possible synonymous sites)
    NP, SP = calculate_np_sp(ref_coding_seq)
    
    # Calculate N (nonsynonymous substitutions) and S (synonymous substitutions)
    N, S = calculate_n_s(ref_coding_seq, target_coding_seq)
    
    # Calculate dN
    try:
        dN = N / NP
    except ZeroDivisionError:
        dN = 0

    # Calculate dN
    try:
        dS = S / SP
    except ZeroDivisionError:
        dS = 0

    # Calculate omega (ω)
    omega = dN / dS

    # Output results
    print(f"NP (possible nonsynonymous sites): {NP}")
    print(f"SP (possible synonymous sites): {SP}")
    print(f"N (nonsynonymous substitutions): {N}")
    print(f"S (synonymous substitutions): {S}")
    print(f"dN (nonsynonymous substitution rate): {dN:.6f}")
    print(f"dS (synonymous substitution rate): {dS:.6f}")
    print(f"ω (dN/dS ratio): {omega:.6f}")
  
    # Add interpretation of omega value here
    if omega > 1:
        print("Interpretation: ω > 1 indicates positive selection.")
    elif omega == 1:
        print("Interpretation: ω = 1 indicates neutral evolution or no selection.")
    elif omega < 1:
        print("Interpretation: ω < 1 indicates purifying selection.")
    else:
        print("Interpretation: Unable to ω value.")

    return omega

# Example usage
if __name__ == '__main__':
    # Define the paths to the input files
    fasta_file = sys.argv[1]
    interval_file = sys.argv[2]
    
    # Run the omega calculation
    calculate_omega(fasta_file, interval_file)
