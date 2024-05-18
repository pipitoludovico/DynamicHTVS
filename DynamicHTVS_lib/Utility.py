import os
from os import listdir


def GetReultFolders(amber) -> list:
    Result_Folders: list = []
    try:
        RESULT_PATH = "post_Docks/" if amber is False else "post_Docks_amber/"
        Result_Folders = [RESULT_PATH + folder for folder in listdir(RESULT_PATH)]
    except:
        print(
            "\nMake sure you have your \"post_Docks\" or \"post_Docks_amber\" folders ready if you want to run the parameterization and dynamics.")
    return Result_Folders


def GetRightSettings() -> None:
    if not os.path.exists('./equilibration'):
        print("Equilibration folder not found. Are you in the right working folder?")
        exit()
