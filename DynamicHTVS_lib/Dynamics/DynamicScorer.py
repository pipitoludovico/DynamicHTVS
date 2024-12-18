from os import getcwd, listdir, path, getpid, cpu_count, system, chdir, makedirs, walk
import sys
import MDAnalysis
import MDAnalysis.analysis.rms
import statistics
from pathlib import Path
from multiprocessing import Pool
from subprocess import Popen
import warnings
import importlib.resources

libdir = str(importlib.resources.files('DynamicHTVS_lib'))
parameter_folder = path.join(libdir, "parameters")

warnings.filterwarnings('ignore')

cpu = int(cpu_count() / 4)
ROOT_ = getcwd()


def writePID():
    mypid = getpid()
    pidFile = open(".mypid_", "w")
    pidFile.write(str(mypid))
    pidFile.close()


def wrap(using_amber=False):
    print("WRAPPING AND FILTERING IN ", getcwd())
    # /scratch/ludovico3/jenny/comparison/l1/vs/crystal/charmm/post_Docks/deprotonated_arachidonic_acid/system/run_1
    if "complex.dcd" not in listdir("../../gbsa"):
        ext = 'dcd' if any(file.endswith('dcd') for file in listdir('./')) else "xtc" if any(file.endswith('xtc') for file in listdir('./')) else ''
        if ext == "":
            print("No trajectory found in", ROOT_)
            return
        trajFile = [path.abspath(file) for file in listdir("./") if file.endswith(ext) and file.startswith("Traj_")][0]
        if using_amber:
            membraneResnames = ('PA', 'ST', 'OL', 'LEO', 'LEN', 'AR', 'DHA', 'PC', 'PE', 'PS', 'PH', 'P2', 'PGR', 'PGS', 'PI', 'CHL')
        else:
            membraneResnames = ("POPC", "PLPC", "PAPE", "POPE", "POPI", "PAPS", "POPA", "SSM", "NSM", "CMH", "CHOL",
                                "DYPC", "YOPC", "PYPE", "YOPE", "POPS", "YOPA", "ERG", "MIPC", "DPPC", "LLPC",
                                "SOPC", "DPPE", "LLPE", "SOPE", "DPPA", "LLPA", "SOPA", "DPPI", "LLPI", "LLPS",
                                "DPPG", "DGDG", "CMH", "SITO", "STIG", "CAMP")
        allMembRes = " ".join(membraneResnames)

        txt = ['package require psfgen\n', 'package require pbctools\n', 'resetpsf\n',
               f'mol new ../structure.psf type psf\n' if not using_amber else 'mol new ../complex.prmtop type parm7\n',
               f'mol addfile ../structure.pdb type pdb\n' if not using_amber else "",
               f'mol addfile {trajFile} type {ext} first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n',
               f'set sel [atomselect top "not (water or ions)"]\n',
               'proc align { rmolid smolid2 seltext } {\n',
               '  set ref [atomselect 0 $seltext frame 0]\n',
               '  set sel [atomselect $smolid2 $seltext]\n',
               '  set all [atomselect $smolid2 all]\n',
               '  set n [molinfo $smolid2 get numframes]\n',
               '  for { set i 1 } { $i < $n } { incr i } {\n',
               '    $sel frame $i\n', '    $all frame $i\n',
               '    $all move [measure fit $sel $ref]\n', '  }\n', '  return\n', '}\n',
               'puts "Unwrapping all"\n'
               'pbc unwrap -all\n'
               'puts "Aligning on frame 0"\n'
               f'align 0 0 "protein and name CA"\n',
               'pbc set {0 0 0} -all\n',
               'animate write dcd ../../gbsa/complex.dcd beg 0 end -1 waitfor all sel $sel\n',
               'quit\n']
        with open("filterTrj.vmd", "w") as vmdscr:
            for line in txt:
                vmdscr.write(line)
        system('vmd -dispdev text -e filterTrj.vmd > filterlog.log 2>&1')


def rmsd(ligName, amber_rmsd):
    data = []
    #/scratch/ludovico3/jenny/comparison/a10/vs/open/post_Docks/deprotonated_oleic_acid/3/system/run_1
    if "RMSDs.dat" not in listdir("../../gbsa"):
        PDB = "../../gbsa/complex.pdb" if not amber_rmsd else "../../gbsa/complex.prmtop"
        XTC = "../../gbsa/complex.dcd"
        u = MDAnalysis.Universe(PDB, XTC)
        ref = MDAnalysis.Universe("../../gbsa/complex.pdb") if not amber_rmsd else u

        R = MDAnalysis.analysis.rms.RMSD(u, ref, select="protein and name CA", groupselections=[f"resname {ligName} and not name H*"])
        R.run(start=0, step=10)
        r_rmsd = R.rmsd.T  # transpose makes it easier for plotting
        data = list(r_rmsd[3])

        with open('../../gbsa/RMSDs.dat', 'w') as rmsdFile:
            for dat in data:
                rmsdFile.write(str(dat) + "\n")
    else:
        with open('../../gbsa/RMSDs.dat', 'r') as rmsdFile:
            for line in rmsdFile:
                data.append(float(line))
    return data


def getPRMTOP(system_=None):
    topFiles = ""
    parFiles = ""
    print(str(parameter_folder))
    for file in listdir(str(parameter_folder)):
        if file.endswith("rtf"):
            topFiles += f"-top {parameter_folder}/{file} "
        if file.endswith('par') or file.endswith('prm'):
            parFiles += f"-param {parameter_folder}/{file} "
        if file.endswith('str'):
            parFiles += f"-str {parameter_folder}/{file} "
    lj_parameter = [CustomParFile for CustomParFile in listdir("../") if CustomParFile.endswith("_LJ.par")][0]
    parmed = open('parmed.inp', 'w')
    txt = ('chamber '
           f'{topFiles}'
           '-top  ../system/new_file_char.top '
           f'{parFiles}'
           f'-param ../{lj_parameter} '
           f'-psf %s.psf '
           f'-crd %s.pdb '
           f'-radii mbondi3\nparmout %s.prmtop') % (
              system_, system_, system_)
    parmed.write(txt)
    parmed.close()
    Popen(f'parmed -i parmed.inp > {system_}.log 2>&1', shell=True).wait()


def write_pbsa_in():
    mmgbsain = open('mmgbsa.in', "w")
    txt = ('&general\n' '\tkeep_files=0, start=1, interval=10\n' '/\n'
           '&gb\n''\tigb=8, saltcon=0.150,\n''/\n')
    mmgbsain.write(txt)
    mmgbsain.close()


def csvTodat() -> list:
    datFile = open("gbsa.dat", "w")
    try:
        data = []
        with open("gbsa.csv", "r") as D:
            for line in D:
                if "DELTA" in line:
                    for deltaLine in D:
                        if ',' in deltaLine:
                            if deltaLine.startswith('Frame'):
                                continue
                            if len(deltaLine.split(',')) == 8:
                                data.append(float(deltaLine.split(',')[-1].strip()))
                                datFile.write(deltaLine.split(',')[-1].strip() + "\n")
        datFile.close()
        return data
    except FileNotFoundError:
        print("No gbsa.csv found. Check if your GBSA analysis went well.")


def gbsa(_amber):
    GBSAs = []
    chdir('../../gbsa')  # /scratch/ludovico3/jenny/comparison/l1/vs/crystal/post_Docks/deprotonated_palmitic_acid/3/gbsa
    try:
        if 'gbsa.csv' not in listdir("./"):
            if _amber is False:
                if any(not path.exists(file) for file in ('complex.prmtop', "receptor.prmtop", 'ligand.prmtop')):
                    systems = ['complex', 'receptor', 'ligand']
                    for i in systems:
                        getPRMTOP(system_=i)
            write_pbsa_in()
            amberPATH = "$AMBERHOME/bin/MMPBSA.py"
            command = f'{amberPATH} -i mmgbsa.in -o results_mmgbsa.dat -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop -y complex.dcd -eo gbsa.csv'
            with open('gbsa.out', 'w') as stdout_file, open('gbsa.err', 'w') as stderr_file:
                Popen(command, shell=True, stdout=stdout_file, stderr=stderr_file).wait()
        if 'gbsa.dat' not in listdir('./'):
            GBSAs = csvTodat()
            return GBSAs
        else:
            if path.getsize("gbsa.dat") != 0:
                with open('gbsa.dat', 'r') as gbsaFile:
                    for gbsaline in gbsaFile:
                        GBSAs.append(float(gbsaline))
                return GBSAs
            else:
                GBSAs = csvTodat()
                return GBSAs
    except Exception as e:
        print("GBSA calculation did not complete. Please inspect your gbsa input files, complex, receptor and ligand files for errors.\n\n")
        print(repr(e))
        print("Using a default array")
        return [1, 1, 1, 1, 1]


def score(RMSDsFull, GBSAsFull):
    if len(RMSDsFull) > 1 and len(GBSAsFull) > 1:
        RMSDs = RMSDsFull[1:]
        GBSAs = GBSAsFull[1:]
        scores = [i / j for i, j in zip(GBSAs, RMSDs)]
        nframes = len(scores)

        scoreSUM = round(sum(scores), 3)
        AvgScore = round(statistics.mean(scores), 3)
        SDScore = round(statistics.pstdev(scores), 3)

        AvgRMSD = round(statistics.mean(RMSDs), 3)
        SDRMSD = round(statistics.pstdev(RMSDs), 3)

        AvgGBSA = round(statistics.mean(GBSAs), 3)
        SDGBSA = round(statistics.pstdev(GBSAs), 3)

        desScore = round((AvgGBSA / AvgRMSD), 3)
        scoresf = open('GBSA_frames.txt', 'w')
        scoresf.write('Frame	Score		GBSA_Energy		RMSD\n')
        for n in range(0, nframes):
            txt = '%s	%s	%s	%s\n' % (str(n), str(scores[n]), str(GBSAs[n]), str(RMSDs[n]))
            scoresf.write(txt)
        scoresf.close()
        with open('DES.txt', 'w') as DESFILE:
            DESFILE.write(
                f"DES: {desScore} SUM: {scoreSUM}, AVG_S: {AvgScore}, STD_S: {SDScore} AVG_GBSA: {AvgGBSA}, SD_GBSA: {SDGBSA} AVG_RMSD: {AvgRMSD}, SD_RMSD: {SDRMSD}")


def getSummary(folder2analize, amber_):

    PATH_ = "DynamicScores"
    if amber_:
        PATH_ += "_amber"
    if path.exists(PATH_):
        system(f'rm {PATH_} -r')
    makedirs(PATH_, exist_ok=True)
    logPath = path.join(ROOT_, PATH_, "DynamicScores.log")
    LOG = open(logPath, 'a')
    for zincFolder, _, files in walk(folder2analize):
        if 'DES.txt' in files:
            scoreslog = Path('%s/DES.txt' % zincFolder)
            print("GETTING WDSF SCORE FROM", scoreslog)  # ok
            if scoreslog.is_file():
                with open(scoreslog, 'r') as f:
                    pathResult = path.abspath(path.join(ROOT_, zincFolder))
                    for line in f:
                        LOG.write(pathResult + " " + line + "\n")
    LOG.close()
    Popen(f'sort {logPath} -k 3n > {PATH_}/sorted_DES.log', shell=True)


def GBSAcalculatorWrapper(folder, mothFolder, gbsa_amber):
    chdir(folder)
    wrap(gbsa_amber)
    RMSDs = rmsd("UNL", gbsa_amber)
    if RMSDs:
        GBSAs = gbsa(gbsa_amber)
        if GBSAs:
            score(RMSDs, GBSAs)
        else:
            print("GBSA failed in", folder)
    else:
        print("RMSD calculation failed in ", folder)
    chdir(mothFolder)


def main(m_amber):
    folder2analize = './post_Docks/' if not m_amber else './post_Docks_amber/'
    print("Folder to analyse ", folder2analize)
    writePID()
    processes = []
    with Pool(processes=cpu) as p:
        for productionFolder, _, files in walk(folder2analize):
            for file in files:
                if file.startswith('Traj_') and file.endswith(('.dcd', '.xtc')):
                    gbsa_process = p.apply_async(GBSAcalculatorWrapper, (productionFolder, ROOT_, amber,))
                    processes.append(gbsa_process)
        for pro in processes:
            try:
                pro.get()
            except Exception as e:
                print(f"GBSA timedout for: {processes.index(pro)}")
                print(e)
                print("Pool will try to complete the remaining processes.")
    getSummary(folder2analize, amber)


if __name__ == '__main__':
    amber = True if sys.argv[1].startswith("T") else False
    main(amber)
