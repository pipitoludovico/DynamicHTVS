import os
import multiprocessing as mp
from subprocess import Popen, DEVNULL, PIPE

from .SmilesCleaner import SmilesCleaner
from .Filter import Filter

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem


class DatabaseDocker:
    def __init__(self, amber, boxsize):
        self.database = None
        self.boxsize = str(boxsize)
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

        except Exception as e:
            print(
                "Database ending in .smi not found. Please provide a database ending in .smi with 'id' and 'smiles' columns")
            print(repr(e))
            exit()
        if not os.path.exists("./equilibration"):
            exit("No folder named 'equilibration' found.")
        else:
            extToCheck = ('.prmtop' if amberCheck else '.psf')
            trjExtensionToCheck = ('.dcd', '.xtc')
            if not any(file.endswith(extToCheck) for file in os.listdir("./equilibration")):
                exit(
                    f"\nRequired files missing in 'equilibration'.\n\nRequired coord and top file extensions: {extToCheck}.")
            if not any(file.endswith(trjExtensionToCheck) for file in os.listdir("./equilibration")):
                exit(f"\n\nMissing required trajecotry with either {trjExtensionToCheck} extension.")

    def GetLastFrameReceptor(self, selection: str, amber: bool) -> [str, str, str]:
        """Get the last frame from your equilibration folder needs to be specifically named this way.
        The function then takes the path of the psf and the pdb that will be used for the docking and the MDs"""
        extensions = ('.parm7', '.psf', '.gro', 'prmtop')
        trajExt = ('.xtc', '.dcd')
        try:
            self.topPath = \
            [os.path.abspath(directory + "/" + file) for directory, _, files in os.walk('./equilibration') for file in
             files if file.endswith(extensions)][0]
            self.trjPath = \
            [os.path.abspath(directory + "/" + file) for directory, _, files in os.walk('./equilibration') for file in
             files if file.endswith(trajExt)][0]
            print("\nUSING: ", self.topPath, " FROM THE EQUILIBRATED SYSTEM AND ", self.trjPath, "AS TRAJECTORY.\n")
        except FileNotFoundError:
            print("No ./equilibration folder found with equilibrated last frame topology and coordinates.")
            exit()

        # check if protein is not ready for docking
        if not os.path.exists('./receptor/system.pdbqt'):
            os.makedirs('receptor', exist_ok=True)
            if self.topPath.endswith('prmtop'):
                ext_ = 'parm7'
            else:
                ext_ = self.topPath.split(".")[-1]

            vmdBuild = [
                f"mol load {ext_} {self.topPath} {self.trjPath.split('.')[-1]} {self.trjPath}\n",
                "package require pbctools\n"
                # 'pbc unwrap\n'
                'set final [atomselect top "not (water or ions)" frame last]\n',
                'set membr [atomselect top "lipid"]\n',
                '$membr set resid [$membr get residue]\n',
                'set all [atomselect top "all" frame last]\n',
                'set protein [atomselect top "protein" frame last]\n',
                '$protein writepdb protein_only.pdb\n'
                '$protein writepsf protein_only.psf\n',
                "$final writepdb system.pdb\n",
                "$final writepsf system.psf\n",
                '$all writepdb last_frame.pdb\n',
                'puts "finished!"\n', "quit\n"]
            with open('last_frame_getter.tcl', 'w') as lastFrGetter:
                for line in vmdBuild:
                    lastFrGetter.write(line)
            # create the system.pdb/psf
            print("\nGetting the last frame...")
            Popen('vmd -dispdev text -e last_frame_getter.tcl > last_frame_getter.log', shell=True).wait()
            # convert to pdbqt with obabel
            print("Converting the receptor with OpenBabel...")
            Popen("obabel -i pdb protein_only.pdb -o pdbqt -O protein.pdbqt -xr", stdout=DEVNULL, stderr=DEVNULL,
                  shell=True).wait()
            print("\n...completed.")
            # remove torsions due obabel faliures
            Popen('grep -Ev "ENDBRANCH|ROOT|ENBRANCH|REMARK|BRANCH|TORSDOF" protein.pdbqt > receptor.pdbqt',
                  shell=True).wait()
            # move the receptor.pdbqt to receptor folder as well as system
            Popen("mv receptor.* receptor; mv system.* receptor;mv last_frame.pdb receptor; mv last_frame_getter.log logs", shell=True).wait()

        # calculating COM once
        COM = ['package require psfgen', 'package require pbctools',
               f'resetpsf\nmol new {self.ROOT}/receptor/system.pdb type pdb',
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
        Popen('rm *.tcl; rm protein.pdbqt protein_only.*', stderr=DEVNULL, stdout=DEVNULL, shell=True).wait()
        return x_coor, y_coor, z_coor

    def MolsFromSmiles(self):
        """Reads a smiles file with a .smi extension and tries to guess the separators.
        The Dataframe is read in chunks by individual processes to speed up the calculation and bypass the GIL.
        The ProcessDF method calls internal functions in parallel as well."""
        os.makedirs('Docking_folder', exist_ok=True)
        try:
            processes = []
            print("Removing stereochemistry and rewriting the SMILES in canonical form...")
            with mp.Pool(processes=self.cpus) as p:
                for df in pd.read_csv(self.database, sep=None, chunksize=1000000, engine='python'):
                    processes.append(p.apply_async(self.ProcessDF, (df,)))
                print("Now dropping duplicates from the dataframe.")
                for proc in processes:
                    proc.get()
            p.close()
            p.join()
        except pd.errors.ParserError:
            print("Couldn't parse your database. Check if your db has the id and smiles columns")

    def ProcessDF(self, df: pd.DataFrame):
        """Canonicalise the SMILES column after dividing the ID from the SMILES. It then removes the stereochemistry
        to avoid redundancies and canonicalise the SMILES. The internal PrepareMols is called in paralle."""
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
        self.PrepareMols(filteredDF, idColName)

    @staticmethod
    def PrepareMols(df: pd.DataFrame, p_idColName: str) -> None:
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

    def DockMolsWrapper(self, d_folder, x_coor, y_coor, z_coor) -> None:
        """Finds the autodock Vina config file that will be used for the parallel dockings.
        No parallelization will occur at this stage as Vina already assigns the CPU resources for an efficient threading.
        Individual folders will be created for each ligand."""
        os.chdir("./Docking_folder/" + d_folder)
        ligand_pdbqt = d_folder + ".pdbqt"
        Popen(f'obabel -i pdb {d_folder}.pdb -o pdbqt -O {ligand_pdbqt} -xn -xh --partialcharge mmff94', stdout=DEVNULL,
              stderr=DEVNULL, shell=True).wait()
        Popen(
            f"qvina2 --center_x {x_coor} --center_y {y_coor} --center_z {z_coor} --size_x {self.boxsize} --size_y {self.boxsize} --size_z {self.boxsize} --receptor {self.ROOT}/receptor/receptor.pdbqt --ligand {ligand_pdbqt} --cpu {self.cpus} --exhaustiveness 32 --num_modes 1 --out {d_folder}_out.pdbqt",
            stdout=DEVNULL, stderr=PIPE,
            shell=True).wait()
        os.chdir(self.ROOT)

    def DockMols(self, selection, amber):
        DockProcesses = []
        x_coor, y_coor, z_coor = self.GetLastFrameReceptor(selection, amber)
        print('Docking coordinates: ', x_coor, y_coor, z_coor)
        with mp.Pool(processes=self.cpus, maxtasksperchild=16) as p:
            for folder in os.listdir('./Docking_folder'):
                DockProcesses.append(p.apply_async(self.DockMolsWrapper, (folder, x_coor, y_coor, z_coor,)))
            for proc in DockProcesses:
                proc.get()
        p.close()
        p.join()

    def GetReceptorPath(self):
        """Returns the path of the receptor."""
        return self.ROOT + "/receptor/system.pdb"
