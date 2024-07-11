import shutil, os

path = "data/prototypes.fasta"

def getprototypeseqs(destinationFolder="RVRefs"):

    # if directory called output not present, create it and copy the file
    if not os.path.exists(destinationFolder):
        os.makedirs(destinationFolder)

    shutil.copyfile(path, f"{destinationFolder}/RVRefs.fasta")

    print(f"The reference sequence have been downloaded to {destinationFolder} directory")
