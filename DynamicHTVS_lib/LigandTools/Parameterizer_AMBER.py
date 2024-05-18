import os

from rdkit import Chem
from multiprocessing import Pool
from subprocess import run, CalledProcessError

import htmd.ui as htmdmodule

cwd = os.getcwd()
cpu = int(os.cpu_count() / 2)


def ParameterizeLigands(p_original_pdb, folder) -> None:
    os.chdir(folder)
    if not os.path.exists('sqm.out') and not os.path.exists('./system'):
        try:
            original_pdb = p_original_pdb.split("_")[0]
            newLigPDBName = str(p_original_pdb).replace('.pdb', '_.pdb')
            newLigMOL2Name = str(p_original_pdb).replace('.pdb', '_.mol2')
            # RDkit to get the formal charge
            mol_for_charge = Chem.MolFromPDBFile(f"../../Docking_folder/{original_pdb}/{original_pdb}.pdb")
            formal_charge = Chem.GetFormalCharge(mol_for_charge)

            ligand = htmdmodule.Molecule(f'{p_original_pdb}')
            ligand.set('segid', 'X')
            ligand.set('resname', 'UNL')
            ligand.write(newLigPDBName)
            del ligand
        except:
            os.chdir(cwd)
            return
        try:
            run(f"antechamber -i {newLigPDBName} -fi pdb -o {newLigMOL2Name} -fo mol2 -s 0 -c bcc -nc {formal_charge} -rn UNL -at gaff2 -pl -1",
                shell=True)
            run(f"parmchk2 -i {newLigMOL2Name} -f mol2 -o UNL.frcmod -s gaff2", shell=True)
        except CalledProcessError:
            print("X" * 200)
            print("antechamber failed with: ", newLigPDBName)
        from DynamicHTVS_lib.LigandTools.Tleap import RunTleap
        RunTleap(newLigMOL2Name, solvate=False, conc=None)
    os.chdir(cwd)


def RunParameterize(ResultFolders) -> None:
    if len(ResultFolders) != 0:
        with Pool(processes=4) as p:
            processes = []
            for folder in ResultFolders:
                for file in os.listdir(folder):
                    if file.endswith('.pdb') and file.startswith(str(folder).split("/")[-1]) and file.endswith(
                            'pose1.pdb'):
                        process = p.apply_async(ParameterizeLigands, (file, folder))
                        processes.append(process)
            for proc in processes:
                proc.get()
        p.close()
        p.join()
    else:
        exit("No post_Docking_amber folder found")
