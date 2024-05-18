import os
from multiprocessing import Pool
from subprocess import DEVNULL
from subprocess import Popen, run, CalledProcessError

cwd = os.getcwd()


def SplitStreamFile(r_file) -> (str, str):
    topName = str(r_file).replace(".pdb", ".top")
    parName = str(r_file).replace(".pdb", ".par")
    top_file = open(topName, 'w')
    par_file = open(parName, 'w')
    count = 0
    intervallo = list()
    streamfileName = f"{str(r_file).replace('.pdb', '')}.str"
    for line in open(streamfileName, 'r'):
        count += 1
        if line[0:3] == 'END':
            intervallo.append(count)
    i = 0
    for line in open(streamfileName, 'r'):
        i += 1
        if i <= int(intervallo[0]):
            top_file.write(line)
        elif (int(intervallo[0]) + 2) < i <= (int(intervallo[1])):
            par_file.write(line)
    top_file.close()
    par_file.close()
    return topName, parName


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


def CalculateRESP2top(cr_file) -> None:
    # splitting .str into par and top
    c_topName, c_parName = SplitStreamFile(cr_file)
    gaussian_output = str(cr_file).replace('.pdb', '.log')
    # running antechamber with the gaussian output to generate the RESP.mol2 file with right charges
    os.system(
        f'/home/scratch/software/amber20/bin/antechamber -i {gaussian_output} -fi gout -o RESP.mol2 -fo mol2 -c resp -eq 2 -s 2 -pf y  -rn UNL')
    os.system('rm esout punch qout QOUT')
    atomNames_charge = []
    try:
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
    except:
        raise Exception("Ligand:", cr_file, " failed parametrization.")


def Parameterize(p_file, folder) -> None:
    os.chdir(folder)
    # entro in bench/post_Docks/ZINC0000232
    os.makedirs('system', exist_ok=True)
    os.makedirs('gbsa', exist_ok=True)

    renamed_file = pdb_rename(p_file)
    gaussian_input = str(p_file).replace('.pdb', '.gau')
    gaussian_check = str(p_file).replace('.pdb', '.chk')
    mol2 = str(p_file).replace(".pdb", ".mol2")
    strFile = str(p_file).replace('.pdb', '.str')
    # building the initial .mol2 file
    run(f'obabel -i pdb {renamed_file} -o mol2 -O {mol2}', stdout=DEVNULL, stderr=DEVNULL, shell=True)
    # using the mol2 file to generate a rough .str file
    try:
        run(f"/home/scratch/software/silcsbio.2022.1/cgenff/cgenff --all-parameters --force-exhaustive {mol2} > {strFile}",
            shell=True)
        with open(strFile, 'r') as cGenFFstreamFile:
            if 'skipped' in cGenFFstreamFile.readlines():
                run(f"/home/scratch/software/silcsbio.2022.1/cgenff/cgenff --all-parameters --force-exhaustive  --bond-order {mol2} > {strFile}",
                    shell=True)
    except:
        raise Exception('Parameterization with SilcBio failed')

    # get the charge for the ligand
    charge = None
    with open(strFile) as streamFile:
        for streamLine in streamFile.readlines():
            if 'RESI' in streamLine:
                charge = int(float(streamLine.split()[2].strip()))
                break
    # getting atom types and coordinates from the renamed file
    atoms, coordinate = [], []
    for line in open(renamed_file):
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
            '%chk=' + gaussian_check + '\n%nprocshared=8\n'
                                       '%mem=8GB\n# HF/LANL2DZ Opt=(Redundant) SCF=Tight Geom=PrintInputOrient pop=(mk,ReadRadii) iop(6/33=2,6/41=10,6/42=17)\n'
                                       '\n<qmtool> simtype="Geometry optimization" </qmtool>\n'
                                       '\n' + f'{charge} {multiplicity}\n')
        for i in coordinate:
            output.write(i)
        for atom in atoms:
            if atom in halogens.keys():
                output.write(f'\n{halogens[atom]}  2.0\n')
    else:
        output.write(
            '%chk=' + gaussian_check + '\n%nprocshared=8\n'
                                       '%mem=8GB\n# HF/6-31G* Opt=(Redundant) SCF=Tight Geom=PrintInputOrient pop=mk iop(6/33=2,6/41=10,6/42=17)\n'
                                       '\n<qmtool> simtype="Geometry optimization" </qmtool>\n'
                                       '\n' + f'{charge} {multiplicity}\n')
        for i in coordinate:
            output.write(i)
        output.write('\n')
    output.close()
    # now we run gaussian to prepare the antechamber input file
    try:
        # sono sempre in bench/post_Docks/ZINC0000232

        run(f'g09 {gaussian_input}', shell=True, check=True)
    except CalledProcessError as e:
        print("Gaussian09 failed with error status:", e)
    else:
        CalculateRESP2top(p_file)
        BuildLJ(p_file)
    finally:

        os.chdir(cwd)
        # torno nel main thread /benchmark/


def BuildLJ(lj_renamed_file) -> None:
    topology = 'new_file_char.top'
    parameters = str(lj_renamed_file).replace(".pdb", ".par")
    toppar = '/home/scratch/MD_utilities/toppar_c36_jul20/par_all36_cgenff.prm'
    atom_type = dict()
    try:
        for line in open(topology):
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
    except:
        raise Exception("Couldn't build LJ parameters for ligand in ", os.getcwd())


def CheckAndRunParameterize() -> None:
    cpu = int(os.cpu_count() / 2)
    for bestPoseFolder in ResultsFolders[:]:
        topology_file = os.path.join(bestPoseFolder, 'new_file_char.top')
        if os.path.exists(topology_file):
            print("Topology already present in in: " + bestPoseFolder + ". Moving to the next ligand.")
            ResultsFolders.remove(bestPoseFolder)
    if len(ResultsFolders) != 0:
        with Pool(processes=4) as p:
            processes = []
            for foldersLeft in ResultsFolders:
                for file in os.listdir(foldersLeft):
                    # trovo in bench/post_Docks/ZINC0000232
                    # folderleft = bench/post_Docks/ZINC0000232
                    if file.endswith('.pdb'):
                        process = p.apply_async(Parameterize, (file, foldersLeft))
                        processes.append(process)
            for proc in processes:
                proc.get()
        p.close()
        p.join()
