from Bio import AlignIO

# Function to compare the lengths of sequences and adjust them to the maximum length by padding with '-'
def compareLengths(seqs):
    max_length = max(len(str(seq)) for seq in seqs)
    adjustedSeqs = [str(seq).ljust(max_length, '-') for seq in seqs]
    return adjustedSeqs

def readFasta(fastaFile):
    alignment = AlignIO.read(fastaFile, "fasta")

    seqList = [record.seq for record in alignment]
    headerList = [record.id for record in alignment]

    seqList = compareLengths(seqList)

    # Combine headers and sequences
    formatted_output = "\n\n".join(f"{header}\n{sequence}" for header, sequence in zip(headerList, seqList))
    return formatted_output
