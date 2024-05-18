from os import makedirs, path
from subprocess import Popen


def FindAndMoveLigands(amber, consider) -> list:
    """Look through the docking results and creates the right folder path"""
    ResultsFolders: list = []
    bestPoses = []
    with open(f'best{consider}.txt', 'r') as bestFile:
        for bestLine in bestFile.readlines():
            idLine = bestLine.split("/")[-1].split("_")[0]
            if idLine not in bestPoses:
                bestPoses.append(idLine)

    for f_pdbPose in bestPoses:
        POST_PATH = "post_Docks"
        if amber:
            POST_PATH += "_amber"
        POST_PATH += "/"

        makedirs(POST_PATH + f_pdbPose, exist_ok=True)
        FULL_RESULT_PATH = path.join(POST_PATH, f_pdbPose)
        ResultsFolders.append(FULL_RESULT_PATH)
        Popen(f'cp ./Docking_folder/{f_pdbPose}/analysis/*_pose1.pdb {FULL_RESULT_PATH}', shell=True).wait()
    return ResultsFolders
