from Bio import SeqIO

# Function 1: Compare lengths among input sequences to pad short sequences with "-" until they are as long as the longest seq
def compare_lengths(seqs):
    # Calculate the maximum length from the list of sequences
    max_length = max(len(seq) for seq in seqs)
    # Adjust each sequence to have the same length as the longest sequence
    adjusted_seqs = [seq.ljust(max_length, '-') for seq in seqs]
    return adjusted_seqs

# Function 2: Function to read sequences from a FASTA file and adjust their lengths
def read_fasta(fasta_file):
    # Read all sequences from the FASTA file using BioPython
    records = list(SeqIO.parse(fasta_file, "fasta"))
    # Extract sequences and headers
    seqs = [str(record.seq).upper() for record in records]
    headers = [record.id for record in records]
    # Adjust all sequences to the length of the longest sequence
    adjusted_seqs = compare_lengths(seqs)
    # Return a dictionary containing the sequences and their corresponding headers
    return {"sequences": adjusted_seqs, "headers": headers}
