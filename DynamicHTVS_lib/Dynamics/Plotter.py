from os import path, getcwd, listdir, makedirs
from multiprocessing.dummy import Pool
import matplotlib.pyplot as plt

import warnings

import MDAnalysis as Mda
import MDAnalysis.analysis.rms

warnings.filterwarnings('ignore')


class DataPlotter:
    def __init__(self, batch):
        self.ROOT: str = getcwd()
        self.batch: int = batch
        self.DynamicPath: str = [path.join(self.ROOT, folder) for folder in listdir(self.ROOT) if "Dynamic" in folder][0]
        self.AMBER = True if "amber" in self.DynamicPath else False
        self.resultTxt = [path.join(self.DynamicPath, textFile) for textFile in listdir(self.DynamicPath) if 'sorted' in textFile][0]
        self.resultsPaths: list = []
        self.sharedDict = {}

    def GetResults(self):
        _ = []
        with open(self.resultTxt, 'r') as resultFile:
            for line in resultFile.readlines():
                if line.split("/")[-3] not in _:
                    _.append(line.split("/")[-3])
                    self.resultsPaths.append(line.split()[0])

    def GetFullRMSD(self, PATH: str):
        PDB = f"{PATH}/complex.pdb" if not self.AMBER else f"{PATH}/complex.prmtop"
        XTC = f"{PATH}/complex.dcd"
        key = PATH.split("/")[-3]
        u = Mda.Universe(PDB, XTC)
        ref = Mda.Universe(f"{PATH}/complex.pdb") if not self.AMBER else u

        R = Mda.analysis.rms.RMSD(u, ref, select="protein and name CA",
                                  groupselections=[f"resname UNL and not name H*"])
        R.run(start=0)
        r_rmsd = R.results.rmsd.T
        self.sharedDict[key] = list(r_rmsd[3])

    def ParallelRMSDs(self):
        self.GetResults()
        processes = []
        with Pool(processes=8) as p:
            for gbsa in self.resultsPaths:
                processes.append(p.apply_async(self.GetFullRMSD, (gbsa,)))
            for proc in processes:
                proc.get()

    def PlotAll(self):
        makedirs("RMSDS", exist_ok=True)
        keys = list(self.sharedDict.keys())
        batches = [keys[i:i + self.batch] for i in range(0, len(keys), self.batch)]
        count = 0
        for batch in batches:
            plt.figure(figsize=(10, 6))
            for k in batch:
                y = self.sharedDict[k]
                X = [frame * 0.2 for frame in range(len(y))]
                if "deprotonated" in k:
                    k = k.split("_")[1]
                plt.style.use('ggplot')  # uses gnuplot style
                plt.plot(X, y, label=k)
            plt.xlabel("Time (ns)")
            plt.ylabel("RMSD (A)")
            plt.title("Comparative Stability Analysis")
            plt.legend()
            plt.savefig(f"RMSDS/RMSD_batch_{count}.png", dpi=300)
