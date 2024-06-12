# rhinotype
rhinotype: A Python Package for the Classification of Rhinoviruses

## Table of contents

1. [Background](#background)
2. [Workflow](#workflow)
3. Installation
4. getprototypeseqs
5. readfasta
6. SNPeek
7. assigntypes
8. pairwisedistances
9. countSNPs
10. plotfrequency
11. plotdistances
12. plottree
13. plotAA
14. citation
15. contributors

## <a id="background"></a>Background 

The Rhinovirus constitutes a significant etiological factor in human respiratory infections, contributing to approximately 50% of the annual incidence of the common cold. This virus is a member of the Enterovirus genus within the Picornaviridae family. The Rhinovirus genus is further subdivided into three species: RV-A, RV-B, and RV-C, which collectively comprise around 169 distinct genotypes. Rhinoviruses are positive-sense, non-enveloped RNA viruses with a genome approximately 7.2 kilobases (kb) in length. The genomic organization includes seven non-structural proteins and four structural proteins—VP1, VP2, VP3, and VP4—integral to the viral replication and host infection mechanisms. Presently, the genotyping of Rhinovirus is conducted manually; however, an R software package referred to as [rhinotypeR](https://github.com/omicscodeathon/rhinotypeR/tree/main), facilitates automated genotyping utilizing the VP4 region. The Python rhinotype package, however, aims to refine genotyping accuracy by incorporating both the VP1 and VP4 regions into the analysis. The VP1 region has demonstrated superior precision in genotypic identification, while the VP4/2 region offers compatibility with a unified amplification protocol, thereby streamlining the genotyping process.

## <a id="workflow"></a>Workflow
