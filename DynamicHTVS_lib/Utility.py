from os import makedirs, path, listdir
from subprocess import Popen


def GetReultFolders(amber) -> list:
    Result_Folders: list = []
    try:
        RESULT_PATH = "post_Docks/" if amber is False else "post_Docks_amber/"
        Result_Folders = [RESULT_PATH + folder for folder in listdir(RESULT_PATH)]
    except FileNotFoundError:
        print(
            "\nMake sure you have your \"post_Docks\" or \"post_Docks_amber\" folders ready if you want to run the parameterization and dynamics.")
    return Result_Folders


def GetRightSettings() -> None:
    if not path.exists('./equilibration'):
        print("Equilibration folder not found. Are you in the right working folder?")
        exit()


def dir_path(string):
    if path.isdir(string):
        return path.abspath(string)
    else:
        raise NotADirectoryError(string)


def file_path(string):
    if path.isfile(string):
        return string
    else:
        raise NotADirectoryError(string)


def check_existence(string):
    if path.exists(string):
        return string
    else:
        raise FileNotFoundError(string)


def dirOrfile(string) -> (str, str):
    response = "dir" if path.isdir(string) else 'smi' if path.isfile(string) and string.endswith(
        'smi') else 'pdb' if path.isfile(string) and string.endswith('pdb') else "smi"
    if response:
        return response, string
    else:
        raise FileNotFoundError('The path you used does not point to any vaild folder, pdb file, or .smi file')


def FindAndMoveLigands(amber, consider) -> list:
    """Look through the docking results and creates the right folder path"""
    POST_PATH = "post_Docks_amber/" if amber else "post_Docks"
    ResultsFolders: list = []
    bestPoses = {}
    # reads the path from the ranked summary
    with open(f'best{consider}.txt', 'r') as bestFile:
        for bestLine in bestFile.readlines():
            bareName = bestLine.split("/")[-1].split("_out")[0]
            poseID = bestLine.split()[0][-1]
            if bareName not in bestPoses:
                bestPoses[bareName] = []
                bestPoses[bareName].append(poseID)
            else:
                bestPoses[bareName].append(poseID)
    for ligandName, poses in bestPoses.items():
        path_ = path.join(POST_PATH, ligandName)
        makedirs(path_, exist_ok=True)
        for pose in poses:
            subpath_ = path.join(f"{POST_PATH}/{ligandName}", pose)
            makedirs(subpath_, exist_ok=True)
            dockingSpecificPoses = path.join("Docking_folder", f"{ligandName}/analysis")
            specificResult = path.join(dockingSpecificPoses, f"{ligandName}_out_pose{pose}.pdb")
            if path_ not in ResultsFolders:
                ResultsFolders.append(path_)
            Popen(f'cp {specificResult} {subpath_}', shell=True).wait()
    return ResultsFolders


def LastFrameWriterCHARMM(topPath: str, trjPath: str) -> None:
    membraneResnames = ["POPC", "POPE", "POPI", "CHL1", "SOPE", "POPS", "SSM"]
    vmdBuild = [
        f"mol load {topPath.split('.')[-1]} {topPath} {trjPath.split('.')[-1]} {trjPath}\n",
        "package require pbctools\n"
        # 'pbc unwrap\n'
        'set final [atomselect top "not (water or ions)" frame last]\n',
        f'set membr [atomselect top "{" ".join(membraneResnames)}"]"\n',
        '$membr set resid [$membr get residue]\n',
        'set all [atomselect top "all" frame last]\n',
        'set protein [atomselect top "protein" frame last]\n',

        '$protein writepdb protein_only.pdb\n'
        '$protein writepsf protein_only.psf\n',

        "$final writepdb forGBSA.pdb\n",
        "$final writepsf forGBSA.psf\n",

        '$all writepdb allAtoms.pdb\n',
        '$all writepsf allAtoms.psf\n',
        'puts "finished!"\n', "quit\n"]
    with open('last_frame_getter.tcl', 'w') as lastFrGetter:
        for line in vmdBuild:
            lastFrGetter.write(line)


def LastFrameWriterAMBER(topPath: str, trjPath: str) -> None:
    trajin_ = [f"parm {topPath}",
               f"trajin {trjPath} lastframe"
               "outtraj allAtoms.pdb"]

    with open('last_frame_getter.in', 'w') as trajinFile:
        for line in trajin_:
            trajinFile.write(line + "\n")


def LastFrameWriterAMBERforGBSA(topPath: str, trjPath: str) -> None:
    trajin_ = [f"parm {topPath}",
               f"trajin {trjPath}",
               "strip :WAT,HOH,TIP3,TIP4P,TIP5P,SPC,SOL",
               "strip :Na+,Cl-,K+,Mg2+,Ca2+,Zn2+,Fe2+,Fe3+,Cu+,Cu2+,Mn2+,Co2+,Ni2+,Br-,I-,Cs+,Rb+,Li+",
               "outtraj forGBSA.pdb pdb"]

    with open('forGBSA.in', 'w') as trajinFile:
        for line in trajin_:
            trajinFile.write(line + "\n")
