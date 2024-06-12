from getprototypeseqs import getprototypeseqs

from readfasta import readFasta

getprototypeseqs(destinationFolder = "RVRefs")

test = readFasta(fastaFile = "RVRefs/RVRefs.fasta")

print(test['sequences'])
