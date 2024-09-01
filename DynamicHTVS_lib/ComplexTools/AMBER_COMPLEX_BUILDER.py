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
    membraneResnames = ('PA', 'ST', 'OL', 'LEO', 'LEN', 'AR', 'DHA', 'PC', 'PE', 'PS', 'PH-', 'P2-', 'PGR', 'PGS', 'PI', 'CHL')
    newLigPDBName = str(p_original_pdb).replace('.pdb', '_.pdb')
    os.chdir(folder)
    os.makedirs('system', exist_ok=True)
    os.makedirs('gbsa', exist_ok=True)
    with open('../../receptor/last_frame.pdb') as complexFile:
        for line in complexFile.readlines():
            if any(memb in line for memb in membraneResnames):
                membraneSystem = True
                break

    if os.path.exists('sqm.out'):
        molfile = [file for file in os.listdir("./") if file.endswith('.mol2')][0]
        if not membraneSystem:
            RECEPTOR_PATH = os.path.abspath('../../receptor/system.pdb')

            new_receptor = htmdmodule.Molecule(RECEPTOR_PATH)
            new_ligand = htmdmodule.Molecule(newLigPDBName)
            complex_ = htmdmodule.Molecule(name='complex')
            complex_.append(new_receptor)
            complex_.append(new_ligand)
            complex_.write('dump.pdb')
            del complex_, new_ligand
            RemoveClashes()
            # preparing the receptor for GBSA => not solvated
            Tleap.RunTleap(RECEPTOR_PATH, solvate=False, conc=None)
            # preparing the ligand for GBSA => not solvated
            Tleap.RunTleap("./complex_noClash.pdb", solvate=False, ionize=False, conc=None, MOL2=molfile)
            numberWaters = run('grep "WAT" solvated.pdb | wc -l', shell=True, capture_output=True, text=True)
            output_string = int(numberWaters.stdout.strip())
            # formula taken from https://computecanada.github.io/molmodsim-amber-md-lesson/12-Adding_Ions/index.html
            ionConc = (0.0028798 * output_string) // 2
            # calculating the waters in the solvated complex
            Tleap.RunTleap('./complex_noClash.pdb', solvate=True, ionize=True, conc=ionConc, MOL2=molfile)
        else:
            # 1 make the ligand + complex from the last frame of the equilibration + ligand pose
            Popen("pdb4amber -i ../../receptor/last_frame.pdb -o last_frame_pdb4amb_connect2.pdb; grep -v 'CONECT' last_frame_pdb4amb_connect2.pdb > last_frame_pdb4amb.pdb", shell=True).wait()
            _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
                 "source leaprc.lipid21", "set default PBRadii mbondi3", "loadamberparams UNL.frcmod",
                 f"UNL = loadmol2 {molfile}", "check UNL", "last_frame = loadpdb last_frame_pdb4amb.pdb",
                 "complex = combine{last_frame UNL}", "savepdb complex complex_initial.pdb", "quit"]
            LocalLeap(_, "initial")
            # we remove clashes from the complex
            RemoveClashes()
            # we add TER as VMD removes them
            run("pdb4amber -i complex_noClash.pdb -o complex_final.pdb", shell=True)
            # assign chains and split into complex, receptor, ligand
            SeparateComponents()
            components = ['complex', 'receptor', 'ligand']
            # we build each component for the GBSA after reassigning the lost TER
            for component in components:
                # component_initial <- comes from the VMD script
                run(f"pdb4amber -i gbsa/{component}_initial.pdb -o gbsa/{component}.pdb", shell=True)
                _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
                     "source leaprc.lipid21", "set default PBRadii mbondi3", "loadamberparams UNL.frcmod",
                     f"UNL = loadmol2 {molfile}", f"{component} = loadpdb gbsa/{component}.pdb",
                     f"saveamberparm {component} ./gbsa/{component}.prmtop ./gbsa/{component}.inpcrd", "quit"]
                LocalLeap(_, component)
            # the production system needs to be the hydrated + noClash and needs the CONECT record removed
            run('grep -v "CONECT" complex_final.pdb > complex_final_NOCONECT.pdb', shell=True)
            _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
                 "source leaprc.lipid21", "set default PBRadii mbondi3", "loadamberparams UNL.frcmod",
                 f"UNL = loadmol2 {molfile}", "check UNL", "complex = loadpdb complex_final_NOCONECT.pdb",
                 'setBox complex "vdw"', "saveamberparm complex ./system/complex.prmtop ./system/complex.inpcrd",
                 'quit']
            LocalLeap(_, "system")
            os.makedirs('logs', exist_ok=True)
            os.makedirs('temp_files', exist_ok=True)
            run('mv *.err logs; mv *.out logs; mv *.log* logs; mv *renum* *nonprot* *sslink temp_files', shell=True)
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
        f"mol load pdb complex_initial.pdb\n",
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


def SeparateComponents():
    membraneResnames = (
        'PA', 'ST', 'OL', 'LEO', 'LEN', 'AR', 'DHA', 'PC', 'PE', 'PS', 'PH', 'P2', 'PGR', 'PGS', 'PI', 'CHL')
    allMembRes = " ".join(membraneResnames)

    vmdSeparate = [
        "package require pbctools\n"
        f"mol load pdb complex_final.pdb\n",
        'set A [atomselect top "protein"]\n',
        '$A set chain P\n',
        f'set B [atomselect top "resname {allMembRes}"]\n',
        '$B set chain M\n',
        'set C [atomselect top "resname UNL"]\n',
        '$C set chain X\n',
        f'set complex [atomselect top "protein or resname {allMembRes} UNL"]\n',
        f'set receptor [atomselect top "protein or resname {allMembRes}"]\n',
        f'set ligand [atomselect top "resname UNL"]\n',
        '$complex writepdb gbsa/complex_initial.pdb\n',
        '$receptor writepdb gbsa/receptor_initial.pdb\n',
        '$ligand writepdb gbsa/ligand_initial.pdb\n',
        'exit\n']
    with open('separate.tcl', 'w') as separator:
        for line in vmdSeparate:
            separator.write(line)
    Popen('vmd -dispdev text -e separate.tcl > separate.log 2>&1', shell=True).wait()


def LocalLeap(_, count):
    with open(f'inleap_{count}', 'w') as inleap:
        for tleapCommand in _:
            inleap.write(tleapCommand + "\n")
    try:
        run(f"tleap -f inleap_{count}", shell=True, stdout=DEVNULL)
    except CalledProcessError:
        print("tleap failed")
        exit()


def Removes():
    pass

# # We first restore TER and then make the receptor.prmtop
# run("pdb4amber -i gbsa/receptor_initial.pdb -o gbsa/receptor.pdb", shell=True)
# _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
#      "source leaprc.lipid21", "set default PBRadii mbondi3",
#      "receptor = loadpdb gbsa/receptor.pdb",
#      "saveamberparm receptor ./gbsa/receptor.prmtop ./gbsa/receptor.inpcrd", "quit"]
# LocalLeap(_, 0)
# # #  3 make the ligand for gbsa
# run("pdb4amber -i gbsa/ligand_initial.pdb -o gbsa/ligand.pdb", shell=True)
# _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
#      "source leaprc.lipid21", "set default PBRadii mbondi3",
#      "loadamberparams UNL.frcmod",
#      f"UNL = loadmol2 {molfile}",
#      "check UNL",
#      "saveamberparm UNL ./gbsa/ligand.prmtop ./gbsa/ligand.inpcrd", "quit"]
# LocalLeap(_, 1)
# # now make the complex for gbsa and simulation
# run("pdb4amber -i gbsa/complex_initial.pdb -o gbsa/complex.pdb", shell=True)
# run('grep -v "CONECT" complex_final.pdb > complex_final_NOCONECT.pdb', shell=True)
# _ = ["source leaprc.protein.ff19SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
#      "source leaprc.lipid21", "set default PBRadii mbondi3",
#      "loadamberparams UNL.frcmod",
#      f"UNL = loadmol2 {molfile}",
#      "check UNL",
#      "dry = loadpdb gbsa/complex.pdb",
#      "saveamberparm dry ./gbsa/complex.prmtop ./gbsa/complex.inpcrd",
#      "complex = loadpdb complex_final_NOCONECT.pdb",
#      'setBox complex "vdw"',
#      "saveamberparm complex ./system/complex.prmtop ./system/complex.inpcrd",
#      'quit']
# LocalLeap(_, 2)
