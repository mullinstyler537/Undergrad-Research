import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO

# Function read_fasta_biopython: 
# Reads sequences from a FASTA file using Biopyton's SeqIO parser.
# Returns a dictionary


def read_fasta_biopython(filename):
    sequences = {}
    for record in SeqIO.parse(filename, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    return sequences


# ---------------------
# Function: gc_content
# ---------------------
# Calculates the GC content (%) of a DNA sequence.
def gc_content(sequence):
    gc = sequence.count('G') + sequence.count('C')
    return (gc / len(sequence)) * 100 if sequence else 0


# ------------------------
# Function: sequence_stats
# ------------------------
# Takes a dictionary of sequences and returns a list of tuples:
# (sequence_id, length, GC_content)
def sequence_stats(sequences):
    stats = []
    for seq_id, seq in sequences.items():
        length = len(seq)
        gc = gc_content(seq)
        stats.append((seq_id, length, gc))
    return stats


# ---------------------------
# Function: find_extremes
# ---------------------------
# Finds the longest and shortest sequences in the dictionary.
# Returns tuples: (seq_id, sequence) for longest and shortest.
def find_extremes(sequences):
    longest = max(sequences.items(), key=lambda x: len(x[1]))
    shortest = min(sequences.items(), key=lambda x: len(x[1]))
    return longest, shortest


# ------------------------
# Function: plot_sequence_lengths
# ------------------------
# Plots a bar chart showing the lengths of each sequence.
def plot_sequence_lengths(stats):
    ids = [item[0] for item in stats]
    lengths = [item[1] for item in stats]
    
    plt.figure(figsize=(10, 5))
    plt.bar(ids, lengths, color='skyblue')
    plt.title("Sequence Lengths")
    plt.xlabel("Sequence ID")
    plt.ylabel("Length (bp)")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


# ------------------------
# Function: plot_gc_content
# ------------------------
# Plots a bar chart showing the GC content (%) of each sequence.
def plot_gc_content(stats):
    ids = [item[0] for item in stats]
    gc_values = [item[2] for item in stats]
    
    plt.figure(figsize=(10, 5))
    plt.bar(ids, gc_values, color='green')
    plt.title("GC Content (%)")
    plt.xlabel("Sequence ID")
    plt.ylabel("GC Content (%)")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()


# ------------------------
# Function: reverse_complement_biopython
# ------------------------
# Returns the reverse complement of a DNA sequence using Biopython.
def reverse_complement_biopython(seq):
    return str(Seq(seq).reverse_complement())


# ------------------------
# Function: find_motif
# ------------------------
# Searches for a DNA motif in each sequence.
# Returns a dictionary {sequence_id: [positions]} where positions are
# zero-based indices where the motif occurs.
def find_motif(sequences, motif):
    motif = motif.upper()
    results = {}
    for seq_id, seq in sequences.items():
        positions = []
        for i in range(len(seq) - len(motif) + 1):
            if seq[i:i+len(motif)] == motif:
                positions.append(i)
        if positions:
            results[seq_id] = positions
    return results


# ------------------------
# Main program
# ------------------------
def main():
    # Ask user for FASTA file name
    print("Paste your FASTA content below. Press Enter, then Ctrl+D (Linux/macOS) or Ctrl+Z (Windows) and Enter to finish:")
    fasta_input = ''
    try:
        while True:
            line = input()
            fasta_input += line + '\n'
    except EOFError:
        pass  # User ends input

    sequences = read_fasta_from_text(fasta_input)

    print(f"\nNumber of sequences: {len(sequences)}\n")
    
    # Calculate and display sequence stats
    stats = sequence_stats(sequences)
    print("Sequence Stats:")
    for seq_id, length, gc in stats:
        print(f"- {seq_id}: Length = {length} bp, GC% = {gc:.2f}")
    
    # Find and display longest and shortest sequences
    longest, shortest = find_extremes(sequences)
    print(f"\nLongest sequence: {longest[0]} ({len(longest[1])} bp)")
    print(f"Shortest sequence: {shortest[0]} ({len(shortest[1])} bp)")
    
    # Plot sequence lengths and GC content
    plot_sequence_lengths(stats)
    plot_gc_content(stats)
    
    # Show reverse complement example (for the first sequence)
    first_seq_id = list(sequences.keys())[0]
    rc = reverse_complement_biopython(sequences[first_seq_id])
    print(f"\nReverse complement of {first_seq_id}:")
    print(rc)
    
    # Search for motif if user wants
    search = input("\nWould you like to search for a motif? (y/n): ").lower()
    if search == 'y':
        motif = input("Enter the DNA motif (e.g., ATG): ").upper()
        motif_results = find_motif(sequences, motif)
        if motif_results:
            print(f"\nMotif '{motif}' found in:")
            for seq_id, positions in motif_results.items():
                print(f"- {seq_id} at positions {positions}")
        else:
            print(f"Motif '{motif}' not found in any sequence.")
    else:
        print("Motif search skipped.")



if __name__=="__main__":
    main()