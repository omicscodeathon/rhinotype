# rhinotype
rhinotype: A Python Package for the Classification of Rhinoviruses

## Table of contents

1. [Background](#background)
2. [Workflow](#workflow)
3. [Installation](#installation)
4. [getprototypeseqs](#getprototypeseqs)
5. [readfasta](#readfasta)
6. [SNPeek](#snpeek)
7. [assigngenotypes](#assigngenotypes)
8. [pairwisedistances](#pairwisedistances)
9. [overallmeandistance](#overallmeandistance)
10. [countSNPs](#countsnps)
11. [plotfrequency](#plotfrequency)
12. [plotdistances](#plotdistances)
13. [plottree](#plottree)
14. [plotAA](#plotaa)
15. [Citation](#citation)
16. [Contributors](#contributors)

## <a id="background"></a>Background 

Among the major causes of human respiratory infections, rhinovirus causes around 50% of all the year-round cases of common cold. This virus belongs to the Enterovirus genus in the Picornaviridae family. RV-A, RV-B, and RV-C are the three species under the genus Rhinovirus, each having about 169 different genotypes. Rhinoviruses are positive-sense, non-enveloped RNA viruses with approximately 7.2 kilobases long genome. The genome codes for seven non-structural proteins and four structural proteins—VP1, VP2, VP3, and VP4—which are part of the genome encoding viral replication and host infection mechanisms. Genotyping of Rhinovirus, until now, is done manually, while an R software package, [rhinotypeR](https://github.com/omicscodeathon/rhinotypeR/tree/main), allows for fully automated genotyping based on the VP4 region. The Python rhinovirus type package brings out VP1 and VP4 regions that aim to improve acuity in genotyping. The VP1 region has been tested to display better acuity in genotypic identification, while the VP4/2 region ensures compatibility under one universal amplification protocol for rapid genotyping.

## <a id="workflow"></a>Workflow

![Workflow](figures/rhinotype%20workflow.svg)

## <a id="installation"></a>Installation

You can install the package using the following command:

```
pip install rhinotype
```

#### Load the package

```
import rhinotype
```

## <a id="getprototypeseqs"></a>getprotypeseqs

Get the prototype sequences of the Rhinovirus genotypes. Making both VP1 and VP4 sequences available. The user should then combine their sequences together with their sequences and align them using MAFFT or any other alignment tool. Then the user should specify where the prototype sequences should be stored, if not specified the getprototypeseqs command will create RvRefs directory and store protoype sequences, by default.

Example:

```
rhinotype.getprototypeseqs()
```

Using own dataset:

```
rhinotype.getprototypeseqs(destination_folder="path to output folder")
```

## <a id="readfasta"></a>readfasta

Reads the alignments/sequences. Compares the input sequences and pads the short sequences with - until they are long as the longest sequence. The sequences are then stored in a dictionary with the sequence name as the key and the sequence as the value.

Example:

```
test = os.path.join(os.path.dirname(__file__), "test.fasta")
rhinotype.readfasta(fasta_file=test)
```

Using own dataset:

```
rhinotype.readfasta(fasta_file="path to fasta file")
```

## <a id="snpeek"></a>SNPeek

SNPeek function visualizes single nuclotide polymorphisms(SNPs) using the users sequences, relative to a specified reference sequence. To specify a reference sequence, the user should move the sequence to the bottom of the alignment. Substitutions are color coded by nucleotide: A = green, T = red, C = blue, G = yellow. The user can specify the output file name, if not specified the output file will be named "SNPeek.png". By default the legend will be False, if the user wants to include the legend, they should set the legend to True.

Example:

```
rhinotype.SNPeek(fasta_file=test)
```
Using own dataset:

```
rhinotype.SNPeek(fasta_file="path to fasta file",  show_legend=False, output_file="output file name")
```

## <a id="assigngenotypes"></a>assigngenotypes

## <a id="pairwisedistances"></a>pairwisedistances

## <a id="overallmeandistance"></a>overallmeandistance

## <a id="countsnps"></a>countSNPs

## <a id="plotfrequency"></a>plotfrequency

## <a id="plotdistances"></a>plotdistances

## <a id="plottree"></a>plottree

## <a id="plotaa"></a>plotAA

## <a id="citation"></a>Citation

## <a id="contributors"></a>Contributors

1. [Ephantus Wambui](https://github.com/Ephantus-Wambui)

2. [Daniel Okoro](https://github.com/danny6200)

3. [Andrew Acheampong](https://github.com/AcheampongAndy)

4. [Parcelli Jepchirchir](https://github.com/Parcelli)

5. [Manase Aloo](https://github.com/manasealoo)
