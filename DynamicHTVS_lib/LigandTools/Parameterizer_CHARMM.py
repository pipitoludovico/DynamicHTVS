from os import listdir, getcwd, path, system, chdir, makedirs
from multiprocessing import Pool
from subprocess import DEVNULL, PIPE
from subprocess import run, CalledProcessError

from rdkit import Chem

cwd = getcwd()


def SplitStreamFile() -> (str, str):
    newTopName, parName = "top.top", "par.par"
    top_file = open(newTopName, 'w')
    par_file = open(parName, 'w')
    count = 0
    intervallo = list()
    initialSTRfile = "strfile.str"
    for line in open(initialSTRfile, 'r'):
        count += 1
        if line[0:3] == 'END':
            intervallo.append(count)
    i = 0
    for line in open(initialSTRfile, 'r'):
        i += 1
        if i <= int(intervallo[0]):
            top_file.write(line)
        elif (int(intervallo[0]) + 2) < i <= (int(intervallo[1])):
            par_file.write(line)
    top_file.close()
    par_file.close()
    return newTopName, parName


def pdb_rename(pdb_input):
    new_line = None
    pdb_output = "renum_" + pdb_input
    output_open = open(pdb_output, 'w')
    n = 0
    for line in open(pdb_input):
        if not (line.startswith('HETATM') or line.startswith('ATOM')):
            continue
        else:
            if ('CL' in line.split()[2].lstrip() or 'BR' in line.split()[2].lstrip()
                    or 'Cl' in line.split()[2].lstrip() or 'Br' in line.split()[2].lstrip()):
                n += 1
                if n < 10:
                    new_line = line[0:12] + '%s%s ' % (line.split()[2][0:2], str(n)) + line[16:-1] + '\n'
                elif 10 <= n <= 99:
                    new_line = line[0:12] + '%s%s' % (line.split()[2][0:2], str(n)) + line[16:-1] + '\n'
                elif n > 99:
                    new_line = line[0:12] + '%s  ' % line.split()[2][0:2] + line[16:-1] + '\n'
                output_open.write(new_line)
            else:
                n += 1
                if n < 10:
                    new_line = line[0:12] + '%s%s  ' % (line.split()[2][0:1], str(n)) + line[16:-1] + '\n'
                elif 10 <= n <= 99:
                    new_line = line[0:12] + '%s%s ' % (line.split()[2][0:1], str(n)) + line[16:-1] + '\n'
                elif n > 99:
                    new_line = line[0:12] + '%s%s' % (line.split()[2][0:1], str(n)) + line[16:-1] + '\n'
                output_open.write(new_line)
    output_open.close()
    return pdb_output


def CalculateRESP2top(pdbForResp_file) -> None:
    # splitting .str into par and top
    c_topName, c_parName = SplitStreamFile()  # top.top par.par
    gaussian_output = "gau_input.log"
    # running antechamber with the gaussian output to generate the RESP.mol2 file with right charges
    if not path.exists("RESP.mol2"):
        system(
            f'/home/scratch/software/amber20/bin/antechamber -i {gaussian_output} -fi gout -o RESP.mol2 -fo mol2 -c resp -eq 2 -s 2 -pf y  -rn UNL > antechamberlog 2>&1')
    if path.exists("QOUT"):
        system('rm esout punch qout QOUT')  # removing a bit of clutter
    if path.exists("RESP.mol2"):
        atomNames_charge = []
        with open('RESP.mol2', 'r') as f:
            for line in f:
                if '@<TRIPOS>BOND' in line:
                    break
                if len(line.split()) == 9:
                    charge = line.split()[8]
                    atomNames_charge.append(charge)
        # writing the new .top file
        new_top_file = open('new_file_char.top', 'w')
        n = 0
        for line in open(c_topName):
            if 'LP' in line:
                continue
            if line[0:4] == 'RESI':
                new_line = line[0:5] + "UNL" + 7 * ' ' + line[15:71] + '\n'
                new_top_file.write(new_line)
            elif line.startswith('ATOM'):
                new_line = line[0:19] + atomNames_charge[n] + '\n'
                new_top_file.write(new_line)
                n += 1
            elif line[0:4] != 'ATOM':
                new_top_file.write(line)
        new_top_file.close()
    else:
        print("Ligand:", pdbForResp_file, " failed parametrization.")


def WriteGaussian(pdbFile, charge) -> None:
    # getting atom types and coordinates from the renamed file
    gaussian_input = f"gau_input.gau"
    gaussian_check = f"gau_chk.chk"
    atoms, coordinate = [], []
    for line in open(pdbFile):
        if line.startswith('HETATM') or line.startswith('ATOM'):
            splitted = line.split()
            new_line = splitted[2] + ' ' + line[30:54].strip() + '\n'
            atoms.append(splitted[2][0])
            coordinate.append(new_line)
    # writing gaussian input file using the atomtypes and coordinates
    output = open(gaussian_input, 'w')
    halogens = {'F': 9, 'Cl': 17, 'Br': 35, 'I': 53, 'At': 85, 'Ts': 117}
    multiplicity = "1"
    if any(halogen in ('F', 'Cl', 'Br', 'I', 'At', 'Ts') for halogen in atoms):
        output.write(
            '%chk=' + gaussian_check + '\n%nprocshared=8\n%mem=8GB\n# HF/LANL2DZ Opt=(Redundant) SCF=Tight Geom=PrintInputOrient pop=(mk,ReadRadii) iop(6/33=2,6/41=10,6/42=17)\n'
                                       '\n<qmtool> simtype="Geometry optimization" </qmtool>\n'
                                       '\n' + f'{charge} {multiplicity}\n')
        for i in coordinate:
            output.write(i)
        for atom in atoms:
            if atom in halogens.keys():
                output.write(f'\n{halogens[atom]}  2.0\n')
    else:
        output.write(
            '%chk=' + gaussian_check + '\n%nprocshared=8\n%mem=8GB\n# HF/6-31G* Opt=(Redundant) SCF=Tight Geom=PrintInputOrient pop=mk iop(6/33=2,6/41=10,6/42=17)\n'
                                       '\n<qmtool> simtype="Geometry optimization" </qmtool>\n'
                                       '\n' + f'{charge} {multiplicity}\n')
        for i in coordinate:
            output.write(i)
        output.write('\n')
    output.close()


def RunSilcBio(mol2, strFile) -> bool:
    try:
        result = run(f"/home/scratch/software/silcsbio.2022.1/cgenff/cgenff --all-parameters --force-exhaustive  {mol2} > {strFile}", shell=True, stdout=DEVNULL, stderr=PIPE, text=True)
        if "skipped" in result.stderr:
            result_second_attempt = run(f"/home/scratch/software/silcsbio.2022.1/cgenff/cgenff --all-parameters --force-exhaustive  --bond-order {mol2} > {strFile}", shell=True, stdout=DEVNULL, stderr=PIPE, text=True)
            if "skipped" in result_second_attempt.stderr:
                print("Silcsbio failed.")
                return False
        else:
            return True
    except FileNotFoundError:
        print("Stream file cration failed with SilcBio.")
        return False


def RenumbererWrapper(initialPDBinPostDockPoses, folder)->None:
    chdir(folder)
    pdb_rename(initialPDBinPostDockPoses)  # renum_whatever.pdb
    chdir(cwd)


def Parameterize(initialPDBinPOSTdocks, folder) -> None:
    chdir(folder)
    molID = str(folder).split("/")[1]
    makedirs('system', exist_ok=True)
    makedirs('gbsa', exist_ok=True)
    renamed_file = pdb_rename(initialPDBinPOSTdocks)  # renum_whatever.pdb
    mol2 = molID + ".mol2"  # zinc23030032.mol2
    strFile = "strfile.str"
    check = True
    # building the initial .mol2 file
    if not path.exists(mol2):
        print("Building the mol2 ", mol2)
        run(f'obabel -i pdb {renamed_file} -o mol2 -O {mol2}', stdout=DEVNULL, stderr=DEVNULL, shell=True)
    # using the mol2 file to generate a rough .str file
    if any(not path.exists(topFile) for topFile in ('strfile.str', 'new_file_char.top')):
        # we make a check exit code bool because silcbio does not return an exit status of 0...
        check = RunSilcBio(mol2, strFile)
    if not check:
        chdir(cwd)
        return
    if not all(path.exists(file) for file in ("par_LJ.par", "new_file_char.top")):
        try:
            print("Building the top and par files.")
            run(f'obabel -i pdb {initialPDBinPOSTdocks} -o mol2 -O {mol2}', stdout=DEVNULL, stderr=DEVNULL, shell=True)
            mol_for_charge = Chem.MolFromMol2File(f"{mol2}")
            charge = Chem.GetFormalCharge(mol_for_charge)
        except Exception as e:
            print("Molecule ", molID, " failed. Turnin off sanitization. Check the results are they might be wrong or unphysical!")
            print(repr(e))
            mol_for_charge = Chem.MolFromPDBFile(f"{initialPDBinPOSTdocks}", sanitize=False)
            charge = Chem.GetFormalCharge(mol_for_charge)
        if not path.exists('gau_input.gau'):
            print(f"CHARGE for {initialPDBinPOSTdocks} = {charge}")
            WriteGaussian(initialPDBinPOSTdocks, charge)  # gau_input.gau, gau_chk.chk
    if not path.exists('gau_input.log'):  # if gaussian was not completed
        try:
            run(f'g09 gau_input.gau', shell=True, check=True)
        except CalledProcessError as e:
            print("Gaussian09 failed with error status:", e, " running antechamber.")
            chdir(cwd)
            return
    if not all(path.exists(file) for file in ("par_LJ.par", "new_file_char.top", 'system/ligand.psf')):
        with open('gau_input.log', 'r') as gauResults:
            for line in gauResults.readlines()[-1:]:
                if len(line.split()) > 2:
                    if "Normal termination of Gaussian" in line:
                        CalculateRESP2top(initialPDBinPOSTdocks)
                        BuildLJ()
                        CreateLigandPsf()
                    else:
                        print('Parameterization with SilcBio and AnteChamber failed.')
    src_dir = path.abspath(".")
    parent_dir = path.abspath("..")

    if all(path.exists(file) for file in ("par_LJ.par", "new_file_char.top")):
        for x in listdir(parent_dir):
            dest_dir = path.abspath(path.join(parent_dir, x))
            if path.isdir(dest_dir) and dest_dir != src_dir:
                run(f"cp -r *.top *.par \"{dest_dir}\"", shell=True)
    chdir(cwd)


def BuildLJ() -> None:
    newTopology = 'new_file_char.top'
    parameters = "par.par"
    toppar = '/home/scratch/MD_utilities/toppar_c36_jul20/par_all36_cgenff.prm'
    atom_type = dict()
    try:
        for line in open(newTopology):
            if line.startswith('ATOM'):
                atom_type[line[12:18]] = atom_type.get(line[12:18], 0) + 1
        i = -1
        for line in open(toppar):
            i += 1
            if line.startswith('NONBONDED'):
                break
        output = parameters.replace('.par', '_LJ.par')
        output_open = open(output, 'w')
        for line in open(parameters):
            if not line.startswith('END'):
                output_open.write(line)
            else:
                break
        output_open.write(
            'NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -\n!cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5\n')
        for line in (open(toppar).readlines())[i:]:
            if line[0:6] in atom_type:
                output_open.write(line)
            if line.startswith('NBFIX'):
                break
        output_open.write('\nEND\n')
        output_open.close()
    except Exception as e:
        print("Couldn't build LJ parameters for ligand in ", getcwd(), e)


def CreateLigandPsf() -> (str, str):
    ex_pdb = [file for file in listdir("./") if file.startswith('renum') and file.endswith('.pdb')]
    renamedPDBinput = str(ex_pdb[0])
    vmdLigand = ['package require psfgen\n', '\n',
                 'topology /home/scratch/MD_utilities/toppar_c36_jul20/top_all36_cgenff.rtf\n',
                 'topology /home/scratch/MD_utilities/toppar_c36_jul20/top_all36_prot.rtf\n',
                 'topology /home/scratch/MD_utilities/toppar_c36_jul20/top_all36_carb.rtf\n',
                 f'topology new_file_char.top\n', '\n', 'segment X {\n', f'pdb {renamedPDBinput}\n', '}\n',
                 '\n', f'coordpdb {renamedPDBinput} X\n', '\n', '\n', 'regenerate angles dihedrals\n', '\n',
                 'guesscoord\n', '\n',
                 f'writepsf ./gbsa/ligand.psf\n',
                 f'writepsf ./system/ligand.psf\n',
                 f'writepdb ./gbsa/ligand.pdb\n',
                 f'writepdb ./system/ligand.pdb\n',
                 'exit\n']
    with open('ligand_psf_maker.tcl', 'w') as ligandMaker:
        for line in vmdLigand:
            ligandMaker.write(line)
    system('vmd -dispdev text -e ligand_psf_maker.tcl 2>&1 > /dev/null')
    if not path.exists("./system/ligand.psf") and not path.exists("./gbsa/ligand.psf"):
        print("\n********** Ligand ", renamedPDBinput, "failed. Skipping molecule. **********\n")


def RunParameterize_CHARMM(ResultsFolders: list) -> None:
    if len(ResultsFolders) != 0:
        with Pool(processes=4) as p:
            processes = []
            for mainLigandFolder in ResultsFolders:
                for poseFolder in sorted(listdir(mainLigandFolder)):
                    subfolder = f"{mainLigandFolder}/{poseFolder}"
                    for file in listdir(subfolder):
                        if file.endswith('.pdb') and not file.startswith('renum'):
                            processes.append(p.apply_async(RenumbererWrapper, (file, subfolder)))
            for proc in processes:
                proc.get()

            processes = []
            for mainLigandFolder in ResultsFolders:
                for poseFolder in sorted(listdir(mainLigandFolder)[0]):
                    subfolder = f"{mainLigandFolder}/{poseFolder}"
                    for file in listdir(subfolder):
                        if file.endswith('.pdb') and not file.startswith('renum'):
                            processes.append(p.apply_async(Parameterize, (file, subfolder)))
            for proc in processes:
                proc.get()
