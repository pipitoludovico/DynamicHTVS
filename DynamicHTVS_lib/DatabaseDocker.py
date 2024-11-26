import os
import multiprocessing as mp
from os import makedirs
from shutil import copy2
from subprocess import Popen, DEVNULL, PIPE

from .SmilesCleaner import SmilesCleaner
from .Filter import Filter
from .Utility import LastFrameWriterCHARMM, LastFrameWriterAMBER, LastFrameWriterAMBERforGBSA

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


class DatabaseDocker:
    def __init__(self, amber, boxsize, ligandType=None, ligandPath=None, poses=1):
        self.database = None
        self.boxsize = str(boxsize)
        self.poses = poses
        self.ligandType: str = ligandType
        self.ligandPath: str = ligandPath
        if self.ligandType == "smi":
            self.CheckNecessaryFiles(amberCheck=amber)
        self.dfList = []
        self.topPath = None
        self.trjPath = None
        self.ROOT = os.getcwd()
        self.cpus = int(os.cpu_count() / 2)

    def CheckNecessaryFiles(self, amberCheck):
        try:
            self.database = [file for file in os.listdir('./') if file.endswith('.smi')]
            if len(self.database) > 1:
                print("More than one .smi file found. Please keep only the .smi file you want to use")
                exit()
            else:
                self.database = self.database[0]
                print("DATABASE USED:", self.database)
        except Exception as e:
            print(e)
            print("Database ending in .smi not found. Please provide a database ending in .smi with 'id' and 'smiles' columns")
            exit()

        if not os.path.exists("./equilibration"):
            exit("No folder named 'equilibration' found.")
        else:
            extToCheck = ('.prmtop' if amberCheck else '.psf')
            trjExtensionToCheck = ('.dcd', '.xtc')
            if not any(file.endswith(extToCheck) for file in os.listdir("./equilibration")):
                exit(f"\nRequired files missing in 'equilibration'.\n\nRequired coord and top file extensions: {extToCheck}.")
            if not any(file.endswith(trjExtensionToCheck) for file in os.listdir("./equilibration")):
                exit(f"\n\nMissing required trajecotry with either {trjExtensionToCheck} extension.")

    def GetLastFrameReceptor(self, selection: str) -> [str, str, str]:
        """Get the last frame from your equilibration folder needs to be specifically named this way.
        The function then takes the path of the psf and the pdb that will be used for the docking and the MDs"""
        extensions = ('.parm7', '.psf', '.gro', 'prmtop')
        trajExt = ('.xtc', '.dcd')
        try:
            self.topPath = [os.path.abspath(directory + "/" + file) for directory, _, files in os.walk('./equilibration') for file in files if file.endswith(extensions)][0]
            self.trjPath = [os.path.abspath(directory + "/" + file) for directory, _, files in os.walk('./equilibration') for file in files if file.endswith(trajExt)][0]
            print("\nUSING: ", self.topPath, " FROM THE EQUILIBRATED SYSTEM AND ", self.trjPath, "AS TRAJECTORY.\n")
        except FileNotFoundError:
            print("No ./equilibration folder found with equilibrated last frame topology and coordinates.")
            exit()
        amberTop, charmTop, gromacsTop = ('parm7', 'prmtop'), 'psf', 'gro'
        FF = 'AMBER' if self.topPath.endswith(amberTop) else 'CHARMM' if self.topPath.endswith(charmTop) else 'GROMACS' if self.topPath.endswith(gromacsTop) else None
        if not FF:
            raise FileNotFoundError("It was not possible to determine the FF style. No topology found.")
        # check if protein is not ready for docking and extract the last frame according to the FF
        os.makedirs('receptor', exist_ok=True)
        if not os.path.exists('./receptor/system.pdbqt'):
            print("\nGetting the last frame...")
            # we need 3 systems:
            if FF == 'CHARMM':
                # writes protein_only.pdb/psf for docking
                # forGBSA.pdb/psf for GBSA
                # and allAtoms.pdb/psf for MD
                LastFrameWriterCHARMM(topPath=self.topPath, trjPath=self.trjPath)
                Popen('vmd -dispdev text -e last_frame_getter.tcl > ./logs/last_frame_getter.log', shell=True).wait()
            if FF == 'AMBER':
                # writes allAtoms.pdb
                LastFrameWriterAMBER(topPath=self.topPath, trjPath=self.trjPath)
                Popen('cpptraj -i last_frame_getter.in > ./logs/last_frame_getter.log;', shell=True, stdout=DEVNULL, stderr=DEVNULL).wait()
                # writes protein_only.pdb
                if os.path.exists('last_frame.pdb'):
                    Popen('pdb4amber -i allAtoms.pdb -o protein_only.pdb -p', shell=True, stdout=DEVNULL, stderr=DEVNULL).wait()
                # this writes forGBSA.pdb
                LastFrameWriterAMBERforGBSA(topPath=self.topPath, trjPath=self.trjPath)

            # convert to pdbqt with obabel
            print("Converting the receptor with OpenBabel for docking...")
            Popen("obabel -i pdb protein_only.pdb -o pdbqt -O protein.pdbqt -xr", stdout=DEVNULL, stderr=DEVNULL, shell=True).wait()
            print("\n...completed.")
            # remove torsions due obabel faliures
            print("\nRemoving torsions and branches info from the pdbqt -> making receptor.pdbqt")
            Popen('grep -Ev "ENDBRANCH|ROOT|ENBRANCH|REMARK|BRANCH|TORSDOF" protein.pdbqt > receptor.pdbqt', shell=True).wait()
            # move the receptor.pdbqt to receptor folder as well as system
            print('Cleaning folders...')
            if FF == 'CHARMM':
                Popen("mv receptor.* receptor; mv *.pdb *.psf receptor;", shell=True).wait()
            if FF == 'AMBER':
                Popen("mv receptor.* receptor; mv *.pdb *.psf receptor; mv *.txt *_sslink *.in receptor;mv last_frame.pdb receptor; mv last_frame_getter.log logs", shell=True).wait()

        # calculating COM once
        COM = ['package require psfgen', 'package require pbctools',
               f'resetpsf\nmol new {self.ROOT}/receptor/allAtoms.pdb type pdb',
               f'set sel [atomselect top "{selection}"]',
               'set geo [measure center $sel]', 'puts "GEO $geo"', 'quit']
        vmdin = Popen(['/usr/local/bin/vmd'], stdin=PIPE, stdout=PIPE, stderr=PIPE, text=True)

        x_coor, y_coor, z_coor = None, None, None
        try:
            for line in COM:
                vmdin.stdin.write(line + '\n')
            vmdin.stdin.close()
            vmdin.wait()
            output_lines = vmdin.stdout.readlines()
            for line in output_lines:
                if "GEO" in line:
                    coords_str = line.strip('\n').split(" ")
                    x_coor, y_coor, z_coor = coords_str[1], coords_str[2], coords_str[3]
                    break
        finally:
            vmdin.stdout.close()
            vmdin.stderr.close()
        # clean ROOT
        Popen('mv *.tcl logs; rm protein.pdbqt protein_only.*', stderr=DEVNULL, stdout=DEVNULL, shell=True).wait()
        return x_coor, y_coor, z_coor

    def MolsFromSmiles(self):
        """Reads a smiles file with a .smi extension and tries to guess the separators.
        The Dataframe is read in chunks by individual processes to speed up the calculation and bypass the GIL.
        The ProcessDF method calls internal functions in parallel as well."""
        try:
            processes = []
            print("Removing stereochemistry and rewriting the SMILES in canonical form...")
            with mp.Pool(processes=self.cpus) as p:
                for df in pd.read_csv(self.database, sep=None, chunksize=1000000, engine='python'):
                    processes.append(p.apply_async(self.WriteCanonicalAndMakePDBS, (df,)))
                print("Now dropping duplicates from the dataframe.")
                for proc in processes:
                    proc.get()
            p.close()
            p.join()
        except Exception as e:
            print("Couldn't parse your database. Check if your db has the id and smiles columns", repr(e))

    def WriteCanonicalAndMakePDBS(self, df: pd.DataFrame):
        """Canonicalise the SMILES column after dividing the ID from the SMILES. It then removes the stereochemistry
        to avoid redundancies and canonicalise the SMILES. The internal PrepareMols is called in parallel."""
        smileCol = [col for col in df.columns if col.lower().startswith('smi')]
        idCol = [col for col in df.columns if "id" in col.lower()]
        initialSmiles = df[smileCol[0]]
        idColName = idCol[0]
        cleaner = SmilesCleaner(initialSmiles)
        canonicalSmiles = cleaner.getCanonicalSmiles()
        df['CanonicalSmiles'] = canonicalSmiles
        try:
            filtered = Filter(df, idColName)
        except IndentationError:
            print(
                "Wrong separator or database format. Try to define a different separator with -se or check your database format ['id', 'smiles']")
            exit()
        filteredDF = filtered.getFiltered()
        self.WriteMolsfromSMI(filteredDF, idColName)

    @staticmethod
    def WriteMolsfromSMI(df: pd.DataFrame, p_idColName: str) -> None:
        """Static method that generates the RDkit::MolObj, add hydrogens and 3D coordinates to ultimately build a pdb file
        that will be used for the docking in parallel. Individual folders will be created later to speed up the process."""
        for index, row in df.iterrows():
            mol = Chem.MolFromSmiles(row['CanonicalSmiles'])
            if mol is not None:
                name = str(row[p_idColName])
                os.makedirs("./Docking_folder/" + name, exist_ok=True)
                mol = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol, randomSeed=42)
                path = f"Docking_folder/{row[p_idColName]}/{row[p_idColName]}.pdb"
                AllChem.MolToPDBFile(mol, path)

    def DockWrapper(self, d_folder, x_coor, y_coor, z_coor) -> None:
        """Finds the autodock Vina config file that will be used for the parallel dockings.
        No parallelization will occur at this stage as Vina already assigns the CPU resources for an efficient threading.
        Individual folders will be created for each ligand."""
        os.chdir("./Docking_folder/" + d_folder)
        ligand_pdbqt = d_folder + ".pdbqt"
        Popen(f'obabel -i pdb {d_folder}.pdb -o pdbqt -O {ligand_pdbqt} -xn -xh --partialcharge mmff94', stdout=DEVNULL, stderr=DEVNULL, shell=True).wait()
        Popen(f"qvina2 --center_x {x_coor} --center_y {y_coor} --center_z {z_coor} --size_x {self.boxsize} --size_y {self.boxsize} --size_z {self.boxsize} --receptor {self.ROOT}/receptor/receptor.pdbqt --ligand {ligand_pdbqt} --exhaustiveness 32 --num_modes {self.poses} --out {d_folder}_out.pdbqt", stdout=DEVNULL, stderr=PIPE, shell=True).wait()
        os.chdir(self.ROOT)

    def DockMols(self, selection):
        print("\nDocking the database...")

        os.makedirs('Docking_folder', exist_ok=True)
        x_coor, y_coor, z_coor = self.GetLastFrameReceptor(selection)
        print('Docking coordinates: ', x_coor, y_coor, z_coor)

        if self.ligandType == "smi":
            self.MolsFromSmiles()  # we make Docking_folder/PDBNAME1/PDBNAME1.pdb for each id in the smiles file
        if self.ligandType == "dir":
            print("\nWARNING: MOLECULES IN THE PDB FORMAT MIGHT NEED TO BE PROTONATED FOR CORRECT CHARGE ASSIGNMENT.\n"
                  "MAKE SURE YOU PROTONATED YOUR LIGANDS RIGHT BEFORE PROCEEDING.\n")
            for file in os.listdir(self.ligandPath):
                newFolderName = f"Docking_folder/{file.split('.')[0]}"
                makedirs(newFolderName, exist_ok=True)
                copy2(os.path.join(self.ligandPath, file), newFolderName)
        if self.ligandType == 'pdb':
            newFolderName = f"Docking_folder/{self.ligandPath.split('.')[0]}"
            makedirs(newFolderName, exist_ok=True)
            copy2(os.path.abspath(self.ligandPath), newFolderName)

        # now we have all the folder/pdb like molecule1/molecule1.pdb...
        DockProcesses = []
        with mp.Pool(processes=self.cpus) as p:
            for folder in os.listdir('./Docking_folder'):
                DockProcesses.append(p.apply_async(self.DockWrapper, (folder, x_coor, y_coor, z_coor,)))
            for proc in DockProcesses:
                proc.get()
        p.close()
        p.join()

    def GetReceptorPath(self):
        """Returns the path of the receptor."""
        return self.ROOT + "/receptor/allAtoms.pdb"
