from subprocess import Popen, DEVNULL
from contextlib import closing
import os
from multiprocessing import Pool


class DockingRanker:
    def __init__(self):
        self.ROOT = os.getcwd()
        self.cpuCount = int(os.cpu_count() / 2)

    def ReadVinaScoresWrapper(self, r_dock_folder) -> None:
        """Read the docking results and splits the poses"""
        os.chdir("./Docking_folder/" + r_dock_folder)
        os.makedirs('analysis', exist_ok=True)
        try:
            docking_result = \
                [file_in_dockFolder for file_in_dockFolder in os.listdir('./') if "_out.pdbqt" in file_in_dockFolder][0]
            Popen(
                f"obabel -i pdbqt {docking_result} -o pdb -O analysis/{docking_result.replace('.pdbqt', '')}_pose.pdb -m",
                shell=True, stderr=DEVNULL, stdout=DEVNULL).wait()
            with open(docking_result, 'r') as docking_result_out:
                modelCount = 1
                for line in docking_result_out.readlines():
                    if "RESULT:" in line:
                        resultSummary = open(f'analysis/{docking_result.replace("_out.pdbqt", f"_{modelCount}.txt")}',
                                             'w')
                        resultSummary.write(
                            self.ROOT + "/Docking_folder/" + str(r_dock_folder) + "/analysis/" + docking_result.replace(
                                '.pdbqt', '_pose') + str(modelCount) + f"\t {line.split()[3]}\t")
                        modelCount += 1
                        resultSummary.close()
        except:
            print("Docking ligand ", r_dock_folder, "failed. Removing the folder.")
            os.chdir(self.ROOT)
            Popen(f'rm -r ./Docking_folder/{r_dock_folder}', shell=True, stdout=DEVNULL).wait()
        os.chdir(self.ROOT)

    def ReadVinaScores(self) -> None:
        with closing(Pool(processes=self.cpuCount)) as p:
            ScoreReadProcesses = []
            for dock_folder in os.listdir('./Docking_folder'):
                ScoreReadProcesses.append(p.apply_async(self.ReadVinaScoresWrapper, (dock_folder,)))
            for proc in ScoreReadProcesses:
                proc.get()
        p.close()
        p.join()

    def GetContactsWrapper(self, folder, g_receptorPath):
        os.chdir(self.ROOT + "/Docking_folder/" + folder + "/analysis")
        poses = [pose for pose in os.listdir('./') if pose.endswith(".pdb") and "out_pose" in pose]
        with open(poses[0], 'r') as poseFile:
            resname = None
            for line in poseFile.readlines():
                if 'ATOM' in line or 'HETATM' in line:
                    resname = line.split()[3]
                    break
        for pose in poses:
            poseNum = pose[-5]
            Popen(f'grep "ATOM" {g_receptorPath} > tempComplex.pdb', shell=True).wait()
            Popen(f'grep "ATOM" {pose} >> tempComplex.pdb', shell=True).wait()
            # adding segid and chain X
            g_complex = open(f'complex.pdb', 'w')
            with open('tempComplex.pdb', 'r') as _:
                for r_line in _.readlines():
                    if resname in r_line:
                        g_complex.write(r_line[0:22] + "X" + r_line[23:74] + "X" + r_line[73:-1] + "\n")
                    else:
                        g_complex.write(r_line)
            # Popen('rm tempComplex.pdb', shell=True).wait()
            g_complex.close()
            package_dir = os.path.dirname(os.path.abspath(__file__))
            tcl_script_path = os.path.join(package_dir, 'VMD', 'getContacts.tcl')
            Popen(f'vmd -dispdev text -e {tcl_script_path}', shell=True, stderr=DEVNULL, stdout=DEVNULL).wait()
            Popen(f'cat contacts.int >> {folder}_{poseNum}.txt', shell=True).wait()
        Popen(f'cat *.txt > {folder}_summary.con', shell=True).wait()
        # Popen(f'rm *.int | rm complex.pdb | rm *.txt', shell=True).wait()
        os.chdir(self.ROOT)

    def GetContacts(self, g_receptorPath: str) -> None:
        g_docking_folders = [g_docking_f for g_docking_f in os.listdir('Docking_folder')]
        for i in range(0, len(g_docking_folders), self.cpuCount):
            subGroup = (g_docking_folders[i:i + self.cpuCount])
            GetContacsProcesses = []
            for folder in subGroup:
                with closing(Pool(processes=self.cpuCount)) as p:
                    GetContacsProcesses.append(p.apply_async(self.GetContactsWrapper, (folder, g_receptorPath,)))
            for proc in GetContacsProcesses:
                proc.get()
            print("Batch ", i, "in progress...")
            GetContacsProcesses.clear()
            p.close()
            p.join()

    def CreateSummaryChart(self):
        if os.path.exists('Summary_chart.txt'):
            Popen("rm Summary_chart.txt", shell=True)
        Popen("touch Summary_chart.txt", shell=True)
        for cs_DockFoler in os.listdir("./Docking_folder"):
            Popen(
                f"cat {self.ROOT}/Docking_folder/{cs_DockFoler}/analysis/{cs_DockFoler}_summary.con >> Summary_chart.txt",
                shell=True).wait()

    def SortSummary(self, consider) -> None:
        for s_file in os.listdir(self.ROOT):
            if s_file == 'Summary_chart.txt':
                Popen('sort Summary_chart.txt -k 2n > sortedSummary.txt', shell=True).wait()
                Popen(f'head -{consider} sortedSummary.txt > best{consider}.txt', shell=True).wait()
                Popen(f'mv Summary_chart.txt logs', shell=True).wait()
                break
