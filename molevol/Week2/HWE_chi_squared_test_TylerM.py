"""
This script tests if a site in a FASTA file at a specified position is in Hardy-Weinberg equilibrium.

It reads a FASTA file containing DNA sequences, calculates allele and genotypic frequencies 
at a specified position, and performs a chi-squared test to determine if the observed genotypic 
frequencies match those expected under Hardy-Weinberg equilibrium.

The script handles heterozygous individuals using IUPAC codes, where each code represents a 
combination of two alleles. For example, 'M' represents 'A'/'C'.


Usage:
    python HWE_chi_squared_test.py <fasta_file> <position>

Example:
    python HWE_chi_squared_test.py seq_data_abbrev.fasta 12

This example will check for Hardy-Weinberg equilibrium at position 12 in the sequences 
provided in the "seq_data_abbrev.fasta" file.
"""

import sys  # Import the sys module to handle command-line arguments
from scipy.stats.distributions import chi2  # Import chi2_contingency from scipy for statistical tests

# IUPAC codes for heterozygous genotypes
IUPAC_CODES = {
    'M': ('A', 'C'),  # A/C
    'R': ('A', 'G'),  # A/G
    'W': ('A', 'T'),  # A/T
    'S': ('C', 'G'),  # C/G
    'Y': ('C', 'T'),  # C/T
    'K': ('G', 'T')   # G/T
}

# Function to parse the FASTA file and extract sequences
def parse_fasta(fasta_file):
    sequences = []  # Initialize a list to store sequences
    with open(fasta_file, 'r') as file:  # Open the FASTA file for reading
        seq = ''  # Initialize an empty string to store a sequence
        for line in file:  # Iterate through each line in the file
            line = line.strip()  # Remove leading/trailing whitespace
            if line.startswith('>'):  # Check if the line is a header (starts with '>')
                if seq:  # If a sequence is already collected, add it to the list
                    sequences.append(seq)
                seq = ''  # Reset the sequence string for the next sequence
            else:
                seq += line  # Append the line to the current sequence
        if seq:  # Add the last sequence if it exists
            sequences.append(seq)
    return sequences  # Return the list of sequences


# Function to calculate allele frequencies at a specified position
def get_allele_frequencies(sequences, position):
    alleles = {}  # Initialize a dictionary to count alleles
    for seq in sequences:  # Iterate through each sequence
        allele = seq[position - 1]  # Get the allele at the specified position
        if allele in IUPAC_CODES:  # Check if the allele is a heterozygous IUPAC code
            for base in IUPAC_CODES[allele]:  # If so, count each base in the IUPAC code
                if base in alleles:
                    alleles[base] += 1
                else:
                    alleles[base] = 1
        else:
            if allele in alleles:  # Check if the allele is already in the dictionary
                alleles[allele] += 2  # Increment the count for this allele
            else:
                alleles[allele] = 2  # Initialize the count for this allele
    total_alleles = sum(alleles.values())  # Calculate the total number of alleles
    obs_allele_frequencies = {allele: count / total_alleles for allele, count in alleles.items()}  # Calculate the frequency for each allele

    # *************************************************************************
    # -------------------------------------------------------------------------
    # Print the observed allele frequencies 
    print("Observed allele frequencies: {}".format(str(obs_allele_frequencies)))
    # -------------------------------------------------------------------------
    # *************************************************************************

    return obs_allele_frequencies  # Return the dictionary of allele frequencies

# Function to calculate observed genotypic frequencies
def get_genotypic_frequencies(sequences, position):
    obs_n_genotypes = {}  # Initialize a dictionary to count genotypes
    for seq in sequences:  # Iterate through each sequence
        allele = seq[position - 1]  # Get the allele at the specified position
        if allele in IUPAC_CODES:  # Check if the allele is a heterozygous IUPAC code
            genotype = ''.join(sorted(IUPAC_CODES[allele]))  # Sort the alleles to create a standard genotype key
        else:
            genotype = allele + allele  # Homozygous genotype
        if genotype in obs_n_genotypes:  # Check if the genotype is already in the dictionary
            obs_n_genotypes[genotype] += 1  # Increment the count for this genotype
        else:
            obs_n_genotypes[genotype] = 1  # Initialize the count for this genotype
    total_genotypes = sum(obs_n_genotypes.values())  # Calculate the total number of genotypes
    obs_genot_frequencies = {genotype: count / total_genotypes for genotype, count in obs_n_genotypes.items()}  # Calculate the frequency for each genotype

    # *************************************************************************
    # -------------------------------------------------------------------------
    # Print the observed genotype numbers (number of individuals) here
    print("Observed genotype numbers: {}".format(str(obs_n_genotypes)))
    # -------------------------------------------------------------------------
    # *************************************************************************

    # *************************************************************************
    # -------------------------------------------------------------------------
    # Print the observed genotype frequencies here
    print("Observed genotype frequencies: {}".format(str(obs_genot_frequencies)))
    # -------------------------------------------------------------------------
    # *************************************************************************

    return obs_n_genotypes, obs_genot_frequencies  # Return the dictionary of genotypes and their frequencies


# Function to calculate expected genotypic frequencies based on allele frequencies
def calculate_expected_genotypic_frequencies(allele_freqs, n):
    # Extract the present alleles and their frequencies
    # n = sample size
    alleles = list(allele_freqs.keys())

    # Initialize a dictionary to store expected genotypic frequencies
    exp_genot_frequencies = {}
    exp_n_genotypes = {}
    # Iterate through each allele pair to calculate expected genotypic frequencies
    for i in range(len(alleles)):
        for j in range(i, len(alleles)):
            allele1 = alleles[i]
            allele2 = alleles[j]
            genotype_list = [alleles[i], alleles[j]]
            genotype = ''.join(sorted(genotype_list))
            if i == j:
                # *************************************************************************
                # -------------------------------------------------------------------------
                # THERE IS A BUG HERE THAT YOU NEED TO FIX
                exp_genot_frequencies[genotype] = allele_freqs[allele1] ** 2  # Expected frequency for the homozygous genotype (e.g., AA)
                # -------------------------------------------------------------------------
                # *************************************************************************

            else:
                # *************************************************************************
                # -------------------------------------------------------------------------
                # THERE IS A BUG HERE THAT YOU NEED TO FIX
               exp_genot_frequencies[genotype] = 2 * allele_freqs[allele1] * allele_freqs[allele2]  # Expected frequency for the heterozygous genotype (e.g., AC)
                # -------------------------------------------------------------------------
                # *************************************************************************

            # *************************************************************************
            # -------------------------------------------------------------------------
            exp_n_genotypes[genotype] = exp_genot_frequencies[genotype] * n
            # -------------------------------------------------------------------------
            # *************************************************************************

    # *************************************************************************
    # -------------------------------------------------------------------------
    # Print the expected genotype numbers (number of individuals) here
    print("Expected genotype numbers: {}".format(str(exp_n_genotypes)))
    # -------------------------------------------------------------------------
    # *************************************************************************

    # *************************************************************************
    # -------------------------------------------------------------------------
    # Print the expected genotype numbers (number of individuals) here
    print("Expected genotype frequencies: {}".format(str(exp_genot_frequencies)))
    # -------------------------------------------------------------------------
    # *************************************************************************

    return exp_n_genotypes, exp_genot_frequencies  # Return the dictionary of expected genotypic frequencies

# -------------------------------------------------------------------------

# Function to perform chi-square test to compare observed and expected frequencies
def chi_square_test(observed, expected):
    genotypes = list(observed)
    total = sum(list(observed.values()))  # Calculate the total number of observed genotypes
    
    # Calculate chi-squared statistic
    chi2_stat = 0
    for genotype in genotypes:
        # *************************************************************************
        # -------------------------------------------------------------------------
        # THERE IS A BUG HERE THAT YOU NEED TO FIX
        chi2_stat += ((observed[genotype] - expected[genotype]) ** 2) / expected[genotype]
        # -------------------------------------------------------------------------
        # *************************************************************************

    # *************************************************************************
    # -------------------------------------------------------------------------
    # Print the chi-squared value
    print("Chi-square value: {}".format(chi2_stat))
    # -------------------------------------------------------------------------
    # *************************************************************************


    return chi2_stat  # Return the chi-squared statistic

# -------------------------------------------------------------------------

# Main function to tie everything together
# *** YOU NEED TO PRINT THE P-VALUE *** #
def main(fasta_file, position):
    sequences = parse_fasta(fasta_file)  # Parse the FASTA file to get sequences
    n_seq = len(sequences)
    allele_frequencies = get_allele_frequencies(sequences, position)  # Get allele frequencies at the specified position
    observed_genotypes, observed_genotype_frequencies = get_genotypic_frequencies(sequences, position)  # Get observed genotypic frequencies
    expected_genotypes, expected_genotype_frequencies = calculate_expected_genotypic_frequencies(allele_frequencies, n_seq)  # Calculate expected genotypic frequencies

    chi2_stat = chi_square_test(observed_genotypes, expected_genotypes)  # Calculate the chi-squared value

    # Calculate p-value from chi-squared test (with 1 degree of freedom for HWE test)
    p_value = chi2.sf(chi2_stat, 1)


    # *************************************************************************
    # -------------------------------------------------------------------------
    # Print the p-value value
    print("P-value from chi-squared test: {:.4f}".format(p_value))
    # -------------------------------------------------------------------------
    # *************************************************************************


    # *************************************************************************
    # -------------------------------------------------------------------------
    # Print the interpretation of the p-value
    if p_value < 0.05:
        print("The p-value is less than 0.05: Reject the null hypothesis. The population is likely NOT in Hardy-Weinberg Equilibrium.")
    else:
        print("The p-value is greater than or equal to 0.05: Fail to reject the null hypothesis. The population is likely in Hardy-Weinberg Equilibrium.")

    # -------------------------------------------------------------------------
    # *************************************************************************

# Script entry point
if __name__ == "__main__":
    if len(sys.argv) != 3:  # Check if the correct number of arguments is provided
        print("Usage: python hwe_test.py <fasta_file> <position>")  # Print usage instructions
    else:
        fasta_file = sys.argv[1]  # Get the FASTA file from the command-line arguments
        position = int(sys.argv[2])  # Get the position as an integer from the command-line arguments
        main(fasta_file, position)  # Call the main function with the provided arguments
