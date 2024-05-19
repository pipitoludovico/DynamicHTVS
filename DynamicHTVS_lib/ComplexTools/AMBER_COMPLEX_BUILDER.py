import os
from subprocess import Popen, DEVNULL, CalledProcessError
from subprocess import run
from multiprocessing import Pool
from os import listdir, getcwd
from DynamicHTVS_lib.LigandTools import Tleap

import htmd.ui as htmdmodule

import logging

logging.getLogger('htmd.ui').setLevel(logging.ERROR)

cwd = getcwd()


def CreateComplex(p_original_pdb, folder) -> None:
    membraneSystem = False
    membraneResnames = (
        'PA', 'ST', 'OL', 'LEO', 'LEN', 'AR', 'DHA', 'PC', 'PE', 'PS', 'PH-', 'P2-', 'PGR', 'PGS', 'PI', 'CHL')
    os.chdir(folder)
    with open('../../receptor/last_frame.pdb') as complexFile:
        for line in complexFile.readlines():
            if any(memb in line for memb in membraneResnames):
                membraneSystem = True
                break

    if os.path.exists('sqm.out'):
        molfile = [file for file in os.listdir("./") if file.endswith('.mol2')][0]
        if not membraneSystem:
            RECEPTOR_PATH = os.path.abspath('../../receptor/system.pdb')
            newLigPDBName = str(p_original_pdb).replace('.pdb', '_.pdb')
            new_receptor = htmdmodule.Molecule(RECEPTOR_PATH)
            new_ligand = htmdmodule.Molecule(newLigPDBName)
            complex_ = htmdmodule.Molecule(name='complex')
            complex_.append(new_receptor)
            complex_.append(new_ligand)
            complex_.write('complex.pdb')
            del complex_, new_ligand
            RemoveClashes()
            # preparing the receptor for GBSA => not solvated
            Tleap.RunTleap(RECEPTOR_PATH, solvate=False, conc=None)
            # preparing the ligand for GBSA => not solvated
            Tleap.RunTleap("./complex.pdb", solvate=False, ionize=False, conc=None, MOL2=molfile)
            numberWaters = run('grep "WAT" solvated.pdb | wc -l', shell=True, capture_output=True, text=True)
            output_string = int(numberWaters.stdout.strip())
            # formula taken from https://computecanada.github.io/molmodsim-amber-md-lesson/12-Adding_Ions/index.html
            ionConc = (0.0028798 * output_string) // 2
            # calculating the waters in the solvated complex
            Tleap.RunTleap('./complex.pdb', solvate=True, ionize=True, conc=ionConc, MOL2=molfile)
        else:
            # 1 clean the receptor + membrane from ions
            Popen('grep -v "[+-]" ../../receptor/system.pdb > receptor.pdb', shell=True).wait()
            # 2 make the receptor for gbsa
            _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
                 "source leaprc.lipid21", "set default PBRadii mbondi3",
                 "receptor = loadpdb receptor.pdb",
                 "saveamberparm receptor ./gbsa/receptor.prmtop ./gbsa/receptor.inpcrd", "quit"]
            LocalLeap(_)
            Popen("cp receptor.pdb ./gbsa/receptor.pdb", shell=True).wait()
            # 3 make the ligand for gbsa
            _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
                 "source leaprc.lipid21", "set default PBRadii mbondi3",
                 "loadamberparams UNL.frcmod",
                 f"UNL = loadmol2 {molfile}",
                 "check UNL",
                 "savepdb UNL ./gbsa/ligand.pdb",
                 "saveamberparm UNL ./gbsa/ligand.prmtop ./gbsa/ligand.inpcrd", "quit"]
            LocalLeap(_)
            # 4 make the ligand + complex for gbsa (no water)
            _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
                 "source leaprc.lipid21", "set default PBRadii mbondi3", "loadamberparams UNL.frcmod",
                 f"UNL = loadmol2 {molfile}",
                 "check UNL",
                 "unsolvated_rec = loadpdb receptor.pdb",
                 "complex = combine{unsolvated_rec UNL}",
                 'setBox complex "vdw"',
                 "savepdb complex ./gbsa/complex.pdb",
                 "saveamberparm complex ./gbsa/complex.prmtop ./gbsa/complex.inpcrd",
                 "quit"]
            LocalLeap(_)
            # MD PREP FILE
            # 5 make the comlex + water + ligand
            Popen("pdb4amber -i ../../receptor/last_frame.pdb -o solvated.pdb", shell=True).wait()
            _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
                 "source leaprc.lipid21", "set default PBRadii mbondi3", "loadamberparams UNL.frcmod",
                 "solvated = loadpdb solvated.pdb",
                 f"docked = loadmol2 {molfile}",
                 "check docked",
                 "complex = combine{solvated docked}",
                 'setBox complex "vdw"',
                 "savepdb complex complex_temp.pdb",
                 "quit"]
            LocalLeap(_)
            # Popen("pdb4amber -i complex_temp.pdb -o complex.pdb", shell=True).wait()
            RemoveClashes()
            Popen("pdb4amber -i complex_noClash.pdb -o complex_final.pdb", shell=True).wait()
            _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
                 "source leaprc.lipid21", "set default PBRadii mbondi3",
                 "loadamberparams UNL.frcmod",
                 f"UNL = loadmol2 {molfile}",
                 "check UNL",
                 "complex = loadpdb complex_final.pdb",
                 'setBox complex "vdw"',
                 "saveamberparm complex ./system/complex.prmtop ./system/complex.inpcrd",
                 'quit']
            LocalLeap(_)
    os.chdir(cwd)


def BuildAMBERsystems(ResultsFolders) -> None:
    if len(ResultsFolders) != 0:
        with Pool(processes=16) as p:
            processes = []
            for foldersLeft in ResultsFolders:
                for file in listdir(foldersLeft):
                    if file.endswith('.pdb') and file.startswith(str(foldersLeft).split("/")[-1]) and file.endswith(
                            'pose1.pdb'):
                        processes.append(p.apply_async(CreateComplex, (file, foldersLeft)))
            for proc in processes:
                proc.get()
        p.close()
        p.join()


def RemoveClashes():
    vmdClash = [
        "package require pbctools\n"
        f"mol load pdb complex_temp.pdb\n",
        'set sel [atomselect top "not (same residue as protein or (resname UNL)) and within 1.3 of resname UNL"]\n',
        'set uniqueChainIDs [lsort -unique [$sel get resid]]\n',
        "if {[llength $uniqueChainIDs] == 0} {lappend uniqueChainIDs 000}\n",
        'set to_keep [atomselect top "all and not resid $uniqueChainIDs"]\n',
        '$to_keep writepdb complex_noClash.pdb\n',
        'exit\n']

    with open('clash_remover.tcl', 'w') as clashRemover:
        for line in vmdClash:
            clashRemover.write(line)
    Popen('vmd -dispdev text -e clash_remover.tcl > clash_remover.log 2>&1', shell=True).wait()


def LocalLeap(_):
    with open('inleap', 'w') as inleap:
        for tleapCommand in _:
            inleap.write(tleapCommand + "\n")
    try:
        run("tleap -f inleap", shell=True, stdout=DEVNULL)
    except CalledProcessError:
        print("tleap failed")
        exit()
