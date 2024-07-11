from Bio import AlignIO

# Function 1: Compare lengths among input sequences to pad short sequences with "-" until they are as long as the longest seq
def compare_lengths(seqs):
    max_length = max(len(seq) for seq in seqs)
    adjusted_seqs = [seq.ljust(max_length, '-') for seq in seqs]
    return adjusted_seqs

# Function 2: Function to read sequences from a FASTA file and adjust their lengths
def read_fasta(fasta_file):
    alignment = AlignIO.read(fasta_file, "fasta")
    seqs = [str(record.seq) for record in alignment]
    headers = [record.id for record in alignment]
    adjusted_seqs = compare_lengths(seqs)
    return {"sequences": adjusted_seqs, "headers": headers}
