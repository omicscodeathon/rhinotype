import shutil, os

def getprototypeseqs(destinationFolder="RVRefs"):
    path = os.path.join(os.path.dirname(__file__), '../data/prototypes.fasta')

    # if directory called output not present, create it and copy the file
    if not os.path.exists(destinationFolder):
        os.makedirs(destinationFolder)

    shutil.copyfile(path, os.path.join(os.path.dirname(__file__), f"../{destinationFolder}/RVRefs.fasta"))
    print(f"The reference sequence have been downloaded to {destinationFolder} directory")
