from os import listdir


def PrepareProtein() -> None:
    if "system.psf" in listdir("./receptor"):
        print("System's psf found!")
        return
    else:
        print("No \"system.psf\" found. Make sure VMD completed before continuing")
        exit(1)
