from os import chdir, getcwd, system, listdir, path, makedirs
from subprocess import Popen
from multiprocessing import Pool

cwd = getcwd()


def CheckComponents():
    if not path.exists("./system") or not path.exists("./gbsa"):
        makedirs("gbsa", exist_ok=True)
        makedirs("system", exist_ok=True)
    if not path.exists("./system/ligand.psf") or not path.exists("./system/ligand.pdb"):
        from DynamicHTVS_lib.LigandTools.Parameterizer_CHARMM import CreateLigandPsf
        CreateLigandPsf()
    if not path.exists("./allAtoms.psf") or not path.exists("./allAtoms.pdb"):
        Popen("cp ../../../receptor/allAtoms.p* .;", shell=True)


def CreateComplex(_workingFolder) -> None:
    chdir(_workingFolder)
    CheckComponents()
    membrane = False
    membraneResnames = ("POPC", "PLPC", "PAPE", "POPE", "POPI", "PAPS", "POPA", "SSM", "NSM", "CMH", "CHOL",
                        "DYPC", "YOPC", "PYPE", "YOPE", "POPS", "YOPA", "ERG", "MIPC", "DPPC", "LLPC",
                        "SOPC", "DPPE", "LLPE", "SOPE", "DPPA", "LLPA", "SOPA", "DPPI", "LLPI", "LLPS",
                        "DPPG", "DGDG", "CMH", "SITO", "STIG", "CAMP")
    with open("../../../receptor/allAtoms.pdb", 'r') as structure:
        for line in structure.readlines():
            if any(resname in line for resname in membraneResnames):
                membrane = True
                break

    try:
        system('cp /home/scratch/MD_utilities/acemd_input_files/combine_parameters.pl  .')
        system(f'./combine_parameters.pl /home/scratch/MD_utilities/toppar_c36_jul20/par_all36_cgenff_empty.prm par_LJ.par > ./system/combined_pars.par')
        makedirs('vmd_logs', exist_ok=True)
        makedirs('vmd_steps_charmm', exist_ok=True)
        vmdMerge = [
            ['package require psfgen\n', 'resetpsf\n',
             f'readpsf ../../../receptor/forGBSA.psf\n', f'readpsf ./system/ligand.psf\n', '\n',
             f'coordpdb ../../../receptor/forGBSA.pdb\n', f'coordpdb ./system/ligand.pdb\n', '\n',
             'writepsf ./system/all.psf\n', 'writepdb ./system/all.pdb\n',
             '\n', 'puts "merging complete!"\n', 'exit\n'],
            [  # this chunk will remove clashes
                'package require psfgen\n', 'resetpsf\n',
                f"mol load psf system/all.psf pdb system/all.pdb\n",
                'set sel [atomselect top "not (same residue as protein or (resname UNL)) and within 1.3 of resname UNL"]\n',
                'set uniqueChainIDs [lsort -unique [$sel get resid]]\n',
                "if {[llength $uniqueChainIDs] == 0} {lappend uniqueChainIDs 000}\n",
                'set to_keep [atomselect top "all and not resid $uniqueChainIDs"]\n',
                '$to_keep writepsf ./gbsa/complex.psf\n', '$to_keep writepdb ./gbsa/complex.pdb\n',
                '$to_keep writepdb ./system/all.pdb\n', '$to_keep writepsf ./system/all.psf\n',
                'exit\n'],
            ['package require psfgen\n', 'package require solvate\n',
             'readpsf ./system/all.psf\n', 'coordpdb ./system/all.pdb\n',
             "solvate ./system/all.psf ./system/all.pdb -b 2.4 -t {} -o ./system/solvated\n".format(
                 "15" if not membrane else "0"),
             'exit\n'],
            ['package require autoionize\n',
             'autoionize -psf ./system/solvated.psf -pdb ./system/solvated.pdb -o ./system/structure -sc 0.154\n',
             'exit\n']]

        for x in range(len(vmdMerge)):
            with open(f'vmd_steps_charmm/vmdFile_{x}.tcl', 'w') as vmdFile:
                for line_ in vmdMerge[x]:
                    vmdFile.write(line_)
            system(f'vmd -dispdev text -e vmd_steps_charmm/vmdFile_{x}.tcl > vmd_logs/vmd_log_{x}.log 2>&1')
        del vmdMerge
        system(f"cp ../../../receptor/forGBSA.pdb ./gbsa/receptor.pdb;cp ../../../receptor/forGBSA.psf ./gbsa/receptor.psf;")
        system(f"rm ./system/all* ./system/solvated.* ./system/ligand.p*")
    except Exception as e:
        print("\nAn error has occurred: ", e,
              "\nNo ligand psf/pbd found. Make sure you run the parameterization module first.")
    chdir(cwd)


def BuildCHARMMsystems(ResultsFolders: list) -> None:
    if len(ResultsFolders) != 0:
        with Pool(processes=8) as p:
            processes = []
            for mainLigandFolder in ResultsFolders:
                for poseFolder in listdir(mainLigandFolder):
                    subfolder = f"{mainLigandFolder}/{poseFolder}"
                    if all(path.exists(path.join(subfolder, file)) for file in ("par_LJ.par", "new_file_char.top")):
                        processes.append(p.apply_async(CreateComplex, (subfolder,)))
            for proc in processes:
                proc.get()
