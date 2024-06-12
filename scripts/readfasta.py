from Bio import AlignIO

# Function to compare the lengths of sequences and adjust them to the maximum length by padding with '-'
def compareLengths(seqs):
    # Find the maximum length of the sequences
    max_length = max(len(str(seq)) for seq in seqs)
    # Adjust sequences to the maximum length by converting Seq objects to strings
    adjustedSeqs = [str(seq).ljust(max_length, '-') for seq in seqs]
    return adjustedSeqs

def readFasta(fastaFile):

    alignment = AlignIO.read(fastaFile, "fasta")

    seqList = [record.seq for record in alignment]
    headerList = [record.id for record in alignment]

    seqList = compareLengths(seqList)

    return {'sequences': seqList, 'headers': headerList}
