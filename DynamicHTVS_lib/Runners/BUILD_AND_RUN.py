import os

import GPUtil
from subprocess import Popen
from multiprocessing import Pool

cwd = os.getcwd()


def GetGPUlist(exclude) -> list:
    """Calls GPUtil to evaluate the avaiable GPUs on the machine"""
    GPUs = GPUtil.getGPUs()
    gpu_ids = []

    for availableGPU in GPUs:
        gpu_ids.append(availableGPU.id)
    if len(gpu_ids) != 0:
        print("Available GPUS: ", gpu_ids)
    return gpu_ids


def CheckFailed(amber) -> None:
    """Check if the complex files were built according to the system preparation sanity builds."""
    # purging folders whose parametrization failed
    post_Docks_PATH = "./post_Docks"
    file_to_check = 'new_file_char.top'
    if amber:
        post_Docks_PATH += "_amber"
        file_to_check = "complex.prmtop"
    for folder, _, check_files in os.walk(post_Docks_PATH):
        if folder.startswith(".post_"):
            continue
        if folder == 'system':
            if file_to_check not in os.listdir(folder):
                os.makedirs('Failed', exist_ok=True)
                with open('failed.txt', 'a') as failed:
                    Popen(f'mv {folder} Failed', shell=True).wait()
                    failed.write(folder + "\n")


def RunnerWrapper(_fol, r_gpu, OPENMM_SCRIPT_PATH, amber, excluded, prt, membraneList) -> None:
    if os.path.isdir(_fol + '/system/'):
        membraneL = f"-membraneRestraints {membraneList}" if membraneList else ""
        os.chdir(_fol + "/system")
        exclusion = "-e " + " ".join(*excluded) + " " if excluded else ""
        if os.path.exists('all.pdb'):
            Popen("rm all.*;rm ligand*;rm solvated*;", shell=True)
        if amber:
            Popen(
                f'python -u {OPENMM_SCRIPT_PATH} {exclusion} -pt 3 -et 1 -prt {prt} {membraneL} -gpu {r_gpu} -in 4 -eq EQ -run RUN > openmm.log 2>&1',
                shell=True).wait()
        else:
            Popen(
                f'python -u {OPENMM_SCRIPT_PATH} {exclusion} -up ./combined_pars.par -ut ./new_file_char.top -pt 3 -et 1 -prt {prt} {membraneL} -gpu {r_gpu} -in 4 -eq EQ -run RUN > openmm.log 2>&1',
                shell=True).wait()
    else:
        with open('failed.txt', 'a') as failed:
            Popen(f'mv {_fol} Failed', shell=True).wait()
            failed.write(_fol + "\n")
    os.chdir(cwd)


def RunOpenMMbatches(workingFolders, SCRIPT_PATH, amber, batch, exclude, prt, membraneList):
    """Using the list of the viable folders, divides this array by the number of available GPUs and a list with size==modulo.
    The number of processes is adjusted on the number of GPUs, unless specified elsewhere."""
    ComplexProcesses = []

    CheckFailed(amber)
    GPUlist = GetGPUlist(exclude)
    grouped_folders = [workingFolders[i:i + batch] for i in range(0, len(workingFolders), batch)]

    idx_GPU = 0
    with Pool(processes=batch) as p:
        for group in grouped_folders:
            for folder_ in group:
                ComplexProcesses.append(p.apply_async(RunnerWrapper, (folder_, idx_GPU, SCRIPT_PATH, amber, exclude, prt, membraneList)))
                idx_GPU = (idx_GPU + 1) % len(GPUlist)
        for complexProc in ComplexProcesses:
            complexProc.get()
        p.close()
        p.join()
