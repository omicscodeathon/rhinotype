import shutil, os

# User input for VP1 or VP4/2 sequence
user_query = input("Do you have VP1 or VP4/2 sequence? (Vp1 or Vp4/2): ")
user_input = user_query.capitalize()

def getprototypeseqs(destinationFolder="RVRefs"):
    # Get prototype sequences according to vp region provided
    def user_get_prototype():
        global path
        if user_input == "Vp1":
            path = os.path.join(os.path.dirname(__file__), '../data/vp1_prototypes.fasta')
        elif user_input == "Vp4/2":
            path = os.path.join(os.path.dirname(__file__), '../data/prototypes.fasta')
        else:
            print("Please provide either VP1 or VP4/2 region")

    # Call user get prototype function
    user_get_prototype()

    # if directory called output not present, create it and copy the file
    if not os.path.exists(destinationFolder):
        os.makedirs(destinationFolder)

    shutil.copyfile(path, os.path.join(os.path.dirname(__file__), f"../{destinationFolder}/RVRefs.fasta"))
    print(f"The reference sequence have been downloaded to {destinationFolder} directory")
