from subprocess import Popen, run, DEVNULL
from multiprocessing import Pool
from os import listdir, getcwd, path, chdir, makedirs
from DynamicHTVS_lib.LigandTools import Tleap

cwd = getcwd()


def CreateComplex(folder) -> None:
    chdir(folder)
    makedirs('system', exist_ok=True)
    makedirs('gbsa', exist_ok=True)
    if path.exists('sqm.out') or path.exists("logs/sqm.out"):
        molfile = [file for file in listdir("./") if file.endswith('.mol2')][0]
        RECEPTOR_PATH = path.abspath('../../../receptor/forGBSA.pdb')
        # prepare the system for GBSA
        Tleap.TleapLigand(mol2path=molfile, name="ligand")  # ligand prmtop
        Tleap.TleapReceptor(recPdbPath=RECEPTOR_PATH, name="receptor")  # receptor prmtop
        Tleap.TleapMakeGBSAComplex(recPdbPath=RECEPTOR_PATH, mol2path=molfile, name='GBSAcomplex')  # complex prmtop

        # prepare now the system for cMD
        # we merge the allAtoms with the ligand
        Tleap.TleapMakeComplexCLASH(recPdbPath="../../../receptor/allAtoms.pdb", mol2path=molfile, name="solvate")
        # then we remove clashes and make a new pdb
        RemoveClashes()  # this takes clashed.pdb through VMD, remove clashes and writes complex_noClash.pdb
        """
        # numberWaters = run('grep "WAT" solvated.pdb | wc -l', shell=True, capture_output=True, text=True)
        # output_string = int(numberWaters.stdout.strip())
        # formula taken from https://computecanada.github.io/molmodsim-amber-md-lesson/12-Adding_Ions/index.html
        # ionConc = (0.0028798 * output_string) // 2
        # calculating the waters in the solvated complex
        """
        run("pdb4amber -i complex_noClash.pdb -o complex_tmp.pdb", shell=True, stdout=DEVNULL, stderr=DEVNULL)
        run('grep -v "CONECT" complex_tmp.pdb > complex_final_NOCONECT.pdb', shell=True)
        Tleap.TleapMakeComplexMD(complexNOCLASHpath="complex_final_NOCONECT.pdb", mol2path=molfile, name='MD')
        # Tleap.TleapIonize(mol2path=molfile, complexPdbPath="solvated_noTER.pdb", conc=ionConc, name="complex")
        makedirs('logs', exist_ok=True)
        makedirs('temp_files', exist_ok=True)
        run('mv *.err logs; mv *.out logs; mv *.log* logs; mv *renum* *nonprot* *sslink temp_files', shell=True)
    chdir(cwd)


def BuildAMBERsystems(ResultsFolders) -> None:
    if len(ResultsFolders) != 0:
        with Pool(processes=8) as p:
            processes = []
            for mainLigandFolder in ResultsFolders:
                for poseFolder in listdir(mainLigandFolder):
                    subfolder = f"{mainLigandFolder}/{poseFolder}"
                    if all(path.exists(path.join(subfolder, file)) for file in ("UNL.frcmod", "sqm.out")):
                        processes.append(p.apply_async(CreateComplex, (subfolder,)))
            for proc in processes:
                proc.get(36000)


def RemoveClashes():
    vmdClash = [
        "package require pbctools\n"
        f"mol load pdb ./clashed.pdb\n",
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
    """Deprecated"""
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


def Removes():
    pass
