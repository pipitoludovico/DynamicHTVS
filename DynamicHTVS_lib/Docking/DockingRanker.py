from subprocess import Popen, DEVNULL
from contextlib import closing
import os
from multiprocessing import Pool
import importlib.resources
import re


class DockingRanker:
    package_dir = str(importlib.resources.files('DynamicHTVS_lib'))

    def __init__(self, restrictions: list):
        self.ROOT: str = os.getcwd()
        self.restrictions: list = restrictions
        if not self.restrictions:
            self.restrictions = ['protein', 'resname UNL']
        self.cpuCount: int = int(os.cpu_count() / 4)

    def ReadVinaScoresWrapper(self, r_dock_folder) -> None:
        """Read the docking results and splits the poses"""
        os.chdir("./Docking_folder/" + r_dock_folder)
        os.makedirs('analysis', exist_ok=True)
        try:
            docking_result = [file_in_dockFolder for file_in_dockFolder in os.listdir('./') if "_out.pdbqt" in file_in_dockFolder][0]
            Popen(f"obabel -i pdbqt {docking_result} -o pdb -O analysis/{docking_result.replace('.pdbqt', '')}_pose.pdb -m", shell=True, stderr=DEVNULL, stdout=DEVNULL).wait()
            with open(docking_result, 'r') as docking_result_out:
                modelCount = 1
                for line in docking_result_out.readlines():
                    if "RESULT:" in line:
                        pathAndScore = self.ROOT + "/Docking_folder/" + str(r_dock_folder) + "/analysis/" + docking_result.replace('.pdbqt', '_pose') + str(modelCount) + f"\t {line.split()[3]}\t"
                        outName = f'analysis/{docking_result.replace("_out.pdbqt", f"_{modelCount}.txt")}'
                        resultSummary = open(outName, 'w')
                        resultSummary.write(pathAndScore)
                        modelCount += 1
                        resultSummary.close()
        except Exception as e:
            print(e)
            print("\nDocking ligand ", r_dock_folder, "failed.")
            os.chdir(self.ROOT)
            # Popen(f'rm -r ./Docking_folder/{r_dock_folder}', shell=True, stdout=DEVNULL).wait()
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
        try:
            for idx, pose in enumerate(poses):
                Popen(f'grep "ATOM" {g_receptorPath} > tempComplex_{idx + 1}.pdb', shell=True).wait()
                Popen(f'grep "ATOM" {pose} >> tempComplex_{idx + 1}.pdb', shell=True).wait()
                #  adding segid and chain X
                nameForContact: str = f"tempComplex_{idx + 1}.pdb"
                tcl_script_path = os.path.join(self.package_dir, 'VMD', 'getContacts.tcl')
                #  We now copy the tcl script for counting the contacts
                Popen(f"cp {tcl_script_path} .", shell=True).wait()
                Popen(f"sed -i 's/PLACEHOLDER_0/{nameForContact}/g' getContacts.tcl", shell=True).wait()
                Popen(f"sed -i 's/PLACEHOLDER_1/{self.restrictions[0]}/g' getContacts.tcl", shell=True).wait()
                Popen(f"sed -i 's/PLACEHOLDER_2/{self.restrictions[1]}/g' getContacts.tcl", shell=True).wait()
                Popen(f"sed -i 's/PLACEHOLDER_3/contacts_{idx + 1}.int/g' getContacts.tcl", shell=True).wait()
                # make contacts.int
                Popen(f'vmd -dispdev text -e getContacts.tcl', shell=True, stderr=DEVNULL, stdout=DEVNULL).wait()

                output_file = f"{folder}_summary.con"
                txt_file = f"{folder}_{idx +1}.txt"
                with open(f"contacts_{idx + 1}.int", 'r') as poseContacts:
                    contact_content = poseContacts.read()
                    numContacts = int(contact_content.split(",")[0])
                    if numContacts > 0:
                        print(f"Found {numContacts} contact(s) in {pose} folder: {folder}")  # deprotonated_arachidonic_acid_out_pose1.pdb, depro_arach_acid
                        Popen(f'cat {txt_file} >> {output_file}', shell=True).wait()
                        Popen(f'cat contacts_{idx + 1}.int >> {output_file}', shell=True).wait()
                    else:
                        Popen(f'rm contacts_{idx + 1}.int', shell=True).wait()
            os.chdir(self.ROOT)
        except Exception as e:
            print("Ligand ", folder, "did not produce results.", e)
        os.chdir(self.ROOT)

    def GetContacts(self, g_receptorPath: str) -> None:
        g_docking_folders = [g_docking_f for g_docking_f in os.listdir('Docking_folder')]

        for i in range(0, len(g_docking_folders), self.cpuCount):
            subGroup = g_docking_folders[i:i + self.cpuCount]
            with closing(Pool(processes=self.cpuCount)) as pool:
                tasks = [pool.apply_async(self.GetContactsWrapper, (folder, g_receptorPath)) for folder in subGroup]
                for task, folder in zip(tasks, subGroup):
                    try:
                        task.get(timeout=600)
                    except Exception as e:
                        print(f"{folder} took too long", e)
                        continue
                pool.close()
                pool.join()

    def CreateSummaryChart(self):
        if os.path.exists('Summary_chart.txt'):
            Popen("rm Summary_chart.txt", shell=True)
        Popen("touch Summary_chart.txt", shell=True)
        for cs_DockFoler in os.listdir("./Docking_folder"):
            Popen(f"cat {self.ROOT}/Docking_folder/{cs_DockFoler}/analysis/{cs_DockFoler}_summary.con >> Summary_chart.txt", shell=True).wait()

    def SortSummary(self, consider) -> None:
        for s_file in os.listdir(self.ROOT):
            if s_file == 'Summary_chart.txt':
                input_file = s_file
                output_file = "output.txt"
                best_scores = {}
                with open(input_file, "r") as infile:
                    for line in infile:
                        parts = line.strip().split("\t")
                        if len(parts) < 2:
                            continue
                        path, score_str = parts[0], parts[1]
                        try:
                            score = float(score_str)
                        except ValueError:
                            continue
                        match = re.search(r'Docking_folder/([^/]+)/', path)
                        if match:
                            ligand = match.group(1)
                            if ligand not in best_scores or score < best_scores[ligand][0]:
                                best_scores[ligand] = (score, line.strip())
                with open(output_file, "w") as outfile:
                    for score, line in best_scores.values():
                        outfile.write(line + "\n")

                Popen('sort output.txt -k 2n > sortedSummary.txt', shell=True).wait()
                Popen(f'head -{consider} sortedSummary.txt > best{consider}.txt', shell=True).wait()
                Popen(f'mv Summary_chart.txt output.txt logs', shell=True).wait()
                break
