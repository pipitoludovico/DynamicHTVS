from os import path, getcwd, chdir, cpu_count, listdir, system

from rdkit import Chem
from multiprocessing import Pool
from subprocess import run, CalledProcessError, DEVNULL

cwd = getcwd()
cpu = int(cpu_count() / 2)


def MakeMol2(pdbPath, folder) -> None:
    chdir(folder)  # post_Docks_amber/ligands/1,2,3... per pose
    try:
        # RDkit to get the formal charge
        molOutName = str(pdbPath).replace('.pdb', '.mol2')
        run(f"sed -i '/^ATOM/!d' {pdbPath}", shell=True)
        if not path.exists(molOutName):
            run(f'obabel -i pdb {pdbPath} -o mol2 -O {molOutName}', stdout=DEVNULL, stderr=DEVNULL, shell=True)
    except Exception as e:
        print("An error occurred during the parameterization of: ", pdbPath, e)
        chdir(cwd)
        return
    chdir(cwd)


def SwapCoordinates(template, destination):
    # we're still in 1
    run(f"cp {destination} {destination}_backup", shell=True)
    with open('temp.mol2', 'w') as mol2:
        with open(template, 'r') as template, open(destination, 'r') as coords:
            for line in template:
                mol2.write(line)
                if "<TRIPOS>ATOM" in line:
                    break
            for line in coords:
                if "<TRIPOS>ATOM" in line:
                    break
            for line_template in template:
                if "@<TRIPOS>SUBSTRUCTURE" in line_template:
                    mol2.write(line_template)
                    for substructure_line in template:
                        mol2.write(substructure_line)
                    break
                line_coords = coords.readline()
                if len(line_coords.split()) == 9:
                    atom_id = line_template[0:19]
                    x = line_coords.split()[2]
                    y = line_coords.split()[3]
                    z = line_coords.split()[4]
                    rest = line_template[50:]
                    mol2.write(f"{atom_id} {x}\t{y}\t{z} {rest}")
                else:
                    mol2.write(line_template)
    run(f"cp temp.mol2 {destination}", shell=True)


def ParameterizeLigands(pdbPath, folder) -> None:
    chdir(folder)  # post_Dock_amber/ligand/1
    molOutName = str(pdbPath).replace('.pdb', '.mol2')
    if not path.exists('sqm.out') and not path.exists('./system'):
        try:
            mol_for_charge = Chem.MolFromMol2File(str(molOutName))
            formal_charge = Chem.GetFormalCharge(mol_for_charge)
        except Exception:
            mol_for_charge = Chem.MolFromMol2File(str(molOutName) + "_backup")
            formal_charge = Chem.GetFormalCharge(mol_for_charge)
        try:
            with open('antechamber.out', 'w') as anteOut, open('antechamber.err', 'w') as anteErr:
                run(f"antechamber -i {pdbPath} -fi pdb -o {molOutName} -fo mol2 -s 0 -c bcc -nc {formal_charge} -rn UNL -at gaff2 -pl -1",
                    shell=True, stdout=anteOut, stderr=anteErr)
                run(f"parmchk2 -i {molOutName} -f mol2 -o UNL.frcmod -s gaff2", shell=True, stdout=anteOut,
                    stderr=anteErr)
        except CalledProcessError:
            print("X" * 200)
            print("antechamber failed with: ", pdbPath)
            chdir(cwd)
    # we copy the coordinates of the different poses but import the atom types and bonds of the parameterized mol2
    # I know... very unelegant. But it works.
    if all(path.exists(file) for file in ("UNL.frcmod", "sqm.out")):
        # we get pose1.mol2 as a template
        templateMOL2 = [path.join(".", mol2) for mol2 in listdir(".") if mol2.endswith(".mol2") and "pose" in mol2][0]
        this = getcwd()
        for x in listdir(path.join(this, "../")):
            full_path = path.join(this, "../", x)
            if x != path.basename(this) and path.isdir(full_path):
                destination = [path.join("../", x, mol2) for mol2 in listdir(full_path) if mol2.endswith(".mol2") and "pose" in mol2][0]
                SwapCoordinates(templateMOL2, destination)
    chdir(cwd)


def RunParameterize(ResultFolders: list) -> None:
    if len(ResultFolders) != 0:
        with Pool(processes=4) as p:
            processes = []
            for mainLigandFolder in ResultFolders:
                for poseFolder in sorted(listdir(mainLigandFolder)):
                    resultPath = path.abspath(path.join(mainLigandFolder, poseFolder))
                    for file in listdir(resultPath):
                        if file.endswith('.pdb') and not file.startswith('renum'):
                            pdbPath = path.join(resultPath, file)
                            process = p.apply_async(MakeMol2, (pdbPath, resultPath))
                            processes.append(process)
            for proc in processes:
                proc.get()

            processes = []
            for mainLigandFolder in ResultFolders:
                for poseFolder in sorted(listdir(mainLigandFolder)[0]):
                    resultPath = path.abspath(path.join(mainLigandFolder, poseFolder))
                    for file in listdir(resultPath):
                        if file.endswith('.pdb') and not file.startswith('renum'):
                            pdbPath = path.join(resultPath, file)
                            process = p.apply_async(ParameterizeLigands, (pdbPath, resultPath))
                            processes.append(process)
            for proc in processes:
                proc.get()
    else:
        exit("No post_Docking_amber folder found")
