from getprototypeseqs import getprototypeseqs

from readfasta import read_fasta

from SNPeek import SNPeek

getprototypeseqs(destinationFolder = "RVRefs")

test = read_fasta(fasta_file = "data/test.fasta")

print(test)

fasta_file = read_fasta(fasta_file="data/test.fasta")

SNPeek(fasta_file)
