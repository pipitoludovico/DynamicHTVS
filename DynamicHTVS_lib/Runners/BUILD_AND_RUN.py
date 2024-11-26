from os import listdir, getcwd, walk, makedirs, path, chdir
from time import sleep
import GPUtil
from subprocess import Popen
from multiprocessing import Pool

cwd = getcwd()


def getGPUids(_excluded):
    GPUs = GPUtil.getGPUs()
    gpu_ids = []
    for availableGPU in GPUs:
        gpu_ids.append(availableGPU.id)
    gpu_ids = list(map(str, gpu_ids))
    if _excluded:
        _excluded = list(map(str, _excluded))
        for ex in _excluded:
            if ex in gpu_ids:
                gpu_ids.remove(ex)
        if len(gpu_ids) != 0:
            print("Available GPUS: ", gpu_ids)
        else:
            raise RuntimeError("No GPUs available for computation.")
    return gpu_ids


def CheckFailed(amber) -> None:
    """Check if the complex files were built according to the system preparation sanity builds."""
    # purging folders whose parametrization failed
    post_Docks_PATH = "./post_Docks"
    file_to_check = 'new_file_char.top'
    if amber:
        post_Docks_PATH += "_amber"
        file_to_check = "complex.prmtop"
    for folder, _, check_files in walk(post_Docks_PATH):
        if folder.startswith(".post_"):
            continue
        if folder == 'system':
            if file_to_check not in listdir(folder):
                makedirs('Failed', exist_ok=True)
                with open('failed.txt', 'a') as failed:
                    Popen(f'mv {folder} Failed', shell=True).wait()
                    failed.write(folder + "\n")


def CheckIfCompleted():
    """ Opens the openmm.log, reads the lines and if it finds the finishind line creates the tagfile"""
    with open('openmm.log', 'r') as openmmLog:
        for line in openmmLog.readlines():
            if "Run finished." in line:
                Popen('touch .dontrestart', shell=True).wait()


def RunnerWrapper(_fol, r_gpu, OPENMM_SCRIPT_PATH, amber, excluded, prt, membraneList) -> None:
    """Runs the openmm pipeline inside the folder, equilibrating and producing a trajectory."""
    chdir(_fol + "/system")
    if path.exists('structure.psf') and path.exists('structure.pdb'):
        membraneL = f"-membraneRestraints {membraneList}" if membraneList else ""
        exclusion = "-e " + " ".join(excluded) + " " if excluded else ""
        if path.exists('all.pdb'):
            Popen("rm all.*;rm ligand*;rm solvated*;", shell=True)
        if not path.exists("new_file_char.top"):
            Popen("cp ../new_file_char.top .", shell=True)
        command = f'python -u {OPENMM_SCRIPT_PATH} {exclusion} -pt 3 -et 1 -prt {prt} {membraneL} -gpu {r_gpu} -in 4 -eq EQ -run RUN > openmm.log 2>&1' if amber else f'python -u {OPENMM_SCRIPT_PATH} {exclusion} -up ./combined_pars.par -ut ./new_file_char.top -pt 3 -et 1 -prt {prt} {membraneL} -gpu {r_gpu} -in 4 -eq EQ -run RUN > openmm.log 2>&1'
        if not path.exists('.dontrestart'):
            print("Cleaning previous unfinished trajectories in ", _fol)
            for previousTraj in listdir('./'):
                if previousTraj.endswith("xtc"):
                    Popen(f"rm *.xtc ./run*/*.xtc", shell=True).wait()
                if previousTraj.endswith("dcd"):
                    Popen(f"rm *.dcd ./run*/*.dcd", shell=True).wait()
            Popen(command, shell=True).wait()
            CheckIfCompleted()
        else:
            print("There is a hidden tagfile called .dontrestart in", _fol,
                  "/system. Remove it if you want to restart your simulation.")
            print("WARNING: if you choose so, you will overwrite your previous results!")
    else:
        with open('failed.txt', 'a') as failed:
            Popen(f'mv {_fol} Failed', shell=True).wait()
            failed.write(_fol + "\n")
        chdir(cwd)
    chdir(cwd)


def RunOpenMMbatches(ligandFolders, SCRIPT_PATH, amber, batch, exclude, prt, membraneList):
    """Using the list of the viable folders, divides this array by the number of available GPUs and a list with size==modulo.
    The number of processes is adjusted on the number of GPUs, unless specified elsewhere."""
    CheckFailed(amber)
    GPUlist = getGPUids(exclude)
    if not batch:
        batch = len(GPUlist)
    ALL = []
    for ligand in ligandFolders:
        for pose in listdir(ligand):
            ALL.append(f"{ligand}/{pose}")
    groupedPostDockMainLigandFolders = [ALL[i:i + batch] for i in range(0, len(ALL), batch)]
    print(groupedPostDockMainLigandFolders)
    print("#" * 200)
    print("Running OpenMM")
    count = 0
    with Pool(processes=batch) as p:
        for poses in groupedPostDockMainLigandFolders:
            ComplexProcesses = []
            for pose in poses:
                idx_GPU = GPUlist[count % len(GPUlist)]
                print(pose, "runs on gpu", idx_GPU)
                ComplexProcesses.append(p.apply_async(RunnerWrapper, (pose, idx_GPU, SCRIPT_PATH, amber, exclude, prt, membraneList)))
                count += 1
            for complexProc in ComplexProcesses:
                try:
                    complexProc.get()
                except Exception as e:
                    print("There is an error in ", repr(e))
    print("#" * 200)
    print("Batches completed")
    print("#" * 200)
