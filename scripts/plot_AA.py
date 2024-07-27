import pandas as pd

def readAA(AAfastaFile):
    # Read all lines from the FASTA file
    with open(AAfastaFile, 'r') as file:
        lines = file.readlines()

    # Initialize lists to store sequences and their headers
    seqList = []
    headerList = []

    # Temporary storage for the current sequence being read
    currentSeq = []

    # Iterate through each line of the FASTA file
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            # If currentSeq is not empty, it means we've finished reading a sequence
            if currentSeq:
                # Join all parts of the sequence into one
                fullSeq = ''.join(currentSeq)
                seqList.append(fullSeq)
            # Reset currentSeq for the next sequence
            currentSeq = []
            # Add the header (without the ">" character) to headerList
            headerList.append(line[1:])
        else:
            # If the line is not a header, it's part of the current sequence
            currentSeq.append(line.upper())

    # Add the last sequence to seqList if it exists
    if currentSeq:
        fullSeq = ''.join(currentSeq)
        seqList.append(fullSeq)

    # Adjust all sequences to the length of the longest sequence
    #seqList = compareLengths(seqList)
    
    # Return a dictionary containing the sequences and their corresponding headers
    return {'sequences': seqList, 'headers': headerList}

def plotAA(AAfastaFile, showLegend=False):
    fastaData = readAA(AAfastaFile)  # Read file
    
    sequences = fastaData['sequences']
    seqNames = fastaData['headers']
    proteinLength = max(len(seq) for seq in sequences)
    
    def compareSequences(seqA, seqB):
        seqAChars = list(seqA)
        seqBChars = list(seqB)
        differences = [i for i in range(len(seqAChars)) if seqAChars[i] != seqBChars[i]]
        subsType = [seqBChars[i] for i in differences]
        return pd.DataFrame({'position': differences, 'subsType': subsType})
    
    colorMap = {
        'R': 'red', 'H': 'red', 'K': 'red',      # Positively charged amino acid 
        'D': 'blue', 'E': 'blue',               # Negatively charged amino acid 
        'S': 'green', 'T': 'green', 'N': 'green', 'Q': 'green', # Polar amino acid 
        'A': 'yellow', 'V': 'yellow', 'I': 'yellow', 'L': 'yellow', 'M': 'yellow', 'F': 'yellow', 
        'W': 'yellow', 'P': 'yellow', 'G': 'yellow', 'Y': 'yellow', 'C': 'yellow'  # Nonpolar amino acid 
    }
    
    diffList = []
    for i in range(1, len(sequences)):
        diff = compareSequences(sequences[0], sequences[i])
        diff['color'] = diff['subsType'].map(colorMap).fillna('gray')
        diffList.append(diff)
    
    # Plotting
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(10, 6))
    plt.xlim(1, proteinLength)
    plt.ylim(0.5, len(sequences))
    plt.xlabel(f"Protein Position of {seqNames[-1]}, acting as reference")
    plt.yticks(range(1, len(sequences) + 1), seqNames, rotation=90)
    
    # Plot small vertical bars for each difference using mapped colors
    for i, diff in enumerate(diffList):
        for j in range(len(diff)):
            plt.plot([diff['position'].iloc[j], diff['position'].iloc[j]], 
                     [i + 1 - 0.25, i + 1 + 0.25], color=diff['color'].iloc[j])
    
    if showLegend:  # Add a semi-transparent legend in the top-left corner
        plt.legend(["+ve charged", "-ve charged", "Polar", "Non-polar", "Other"],
                   loc='upper left', bbox_to_anchor=(1, 1), framealpha=0.7)
    
    plt.show()
