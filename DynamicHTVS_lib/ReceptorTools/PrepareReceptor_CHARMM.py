from os import listdir


def PrepareProtein() -> None:
    if "allAtoms.psf" in listdir("./receptor"):
        print("System's psf found!")
        return
    else:
        print("No \"allAtoms.psf\" found. Make sure VMD completed before continuing")
        exit(1)
