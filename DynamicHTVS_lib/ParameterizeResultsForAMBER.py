# import os
# from multiprocessing import Pool
# from subprocess import Popen, run, CalledProcessError, DEVNULL
#
# ResultsFolders: list = []
# cwd = os.getcwd()
# cpu = int(os.cpu_count() / 2)
#
#
# def ParameterizeAndCreateComplex(p_original_pdb, folder, QM=False) -> None:
#     os.chdir(folder)
#     # sono dentro /scratch/ludovico9/galectin/benchmark/post_Docks_amber/ZINC000238950253
#     # se sqm non c'è  e system non è stato creato per la dinamica => no antechamber quindi lo lancio
#     if not os.path.exists('sqm.out') and not os.path.exists('./system'):
#         original_pdb = p_original_pdb.split("_")[0]
#         # newLigPDBName = str(p_original_pdb).replace('.pdb', '_.pdb')
#         # newLigMOL2Name = str(p_original_pdb).replace('.pdb', '_.mol2')
#
#         # mol_for_charge = Chem.MolFromPDBFile(f"../../Docking_folder/{original_pdb}/{original_pdb}.pdb")
#         # formal_charge = Chem.GetFormalCharge(mol_for_charge)
#
#         # new_receptor = Molecule('../../receptor/receptor.pdb')
#
#         # new_ligand = Molecule(newLigPDBName)
#
#         # complex_ = Molecule(name='complex')
#         # complex_.append(new_receptor)
#         # complex_.append(new_ligand)
#         # complex_.write('complex.pdb')
#         # del complex_, new_ligand
#
#         # if QM:
#         #     WriteGaussianInput(newLigPDBName, formal_charge)
#         # else:
#         #     run(f"echo antechamber -i {newLigPDBName} -fi pdb -o {newLigMOL2Name} -fo mol2 -s 0 -c bcc -nc {formal_charge} -rn UNL -at gaff2 -pl 15",
#         #         shell=True)
#         #     run(f"antechamber -i {newLigPDBName} -fi pdb -o {newLigMOL2Name} -fo mol2 -s 0 -c bcc -nc {formal_charge} -rn UNL -at gaff2 -pl 15",
#         #         shell=True)
#         # del mol_for_charge, formal_charge
#
#         # run(f"parmchk2 -i {newLigMOL2Name} -f mol2 -o UNL.frcmod -s gaff2", shell=True)
#
#         # os.makedirs('system', exist_ok=True)
#         # os.makedirs('gbsa', exist_ok=True)
#         # sono sempre in /scratch/ludovico9/galectin/benchmark/post_Docks_amber/ZINC000238950253
#         # _ = ["source leaprc.protein.ff14SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
#         #      "set default PBRadii mbondi2",
#         #      "loadamberparams UNL.frcmod", f"UNL = loadmol2 {newLigMOL2Name}", "check UNL",
#         #      "saveoff UNL UNL.lib", "saveamberparm UNL ./gbsa/ligand.prmtop ./gbsa/ligand.inpcrd",
#         #      'rec = loadpdb ../../receptor/receptor.pdb',
#         #      "saveamberparm rec ./gbsa/receptor.prmtop ./gbsa/receptor.inpcrd",
#         #      "complex = loadpdb complex.pdb",
#         #      "saveamberparm complex ./gbsa/complex.prmtop ./gbsa/complex.inpcrd"
#         #      'setBox complex "vdw"', "solvateBox complex TIP3PBOX 15 0.75",
#         #      "savepdb complex solvated.pdb",
#         #      "quit"]
#         #
#         # with open('inleap', 'w') as inleap:
#         #     for tleapCommand in _:
#         #         inleap.write(tleapCommand + "\n")
#         # try:
#             # creo complex prmtop e inpcrd (NON SOLVATATO) in gbsa per la GBSA,
#             # ligand e receptor prmtop incrd in system e gbsa
#             # run("tleap -f inleap", shell=True, stdout=DEVNULL)
#         # except CalledProcessError:
#         #     print("tleap failed")
#         #     exit()
#         # numberWaters = run('grep "WAT" solvated.pdb | wc -l', shell=True, capture_output=True, text=True)
#         # output_string = int(numberWaters.stdout.strip())
#         # formula taken from https://computecanada.github.io/molmodsim-amber-md-lesson/12-Adding_Ions/index.html
#         # conc = (0.0028798 * output_string) // 2
#
#         # sono in /scratch/ludovico9/galectin/benchmark/post_Docks_amber/ZINC000238950253
#         # _ = ["source leaprc.protein.ff14SB", "source leaprc.gaff2", "source leaprc.water.tip3p",
#         #      "set default PBRadii mbondi2", "loadamberparams UNL.frcmod", f"UNL = loadmol2 {newLigMOL2Name}",
#         #      "check UNL",
#         #      "complex = loadpdb solvated.pdb", "addIons complex Na+ 0",
#         #      "addIons complex Cl- 0", f"addIons complex Na+ {conc}", f"addIons complex Cl- {conc}",
#         #      'setBox complex "vdw"',
#         #      "saveamberparm complex ./system/complex.prmtop ./system/complex.inpcrd",
#         #      "quit"]
#         #
#         # with open('inleap', 'w') as inleap:
#         #     for tleapCommand in _:
#         #         inleap.write(tleapCommand + "\n")
#         # try:
#             # creo il complesso CON LE ACQUE per la dinamica in system
#             # run("tleap -f inleap", shell=True, stdout=DEVNULL)
#         # dentro /scratch/ludovico9/galectin/benchmark/post_Docks_amber/ZINC000238950253/system e gbsa ora ci sono complex e rec e ligand prmtop
#         # solo che in ./gbsa complex non è solvatato o ha ioni
#         # except CalledProcessError:
#         #     print("tleap failed. Check your initial pdb structure and ligand sanity.")
#         #     exit()
#     os.chdir(cwd)
#     # sono nel main thread /benchmark
#
#
# def WriteGaussianInput(newLigMOL2Name, wg_total_charge):
#     gaussian_input = newLigMOL2Name.replace('.mol2', '.gau')
#     gaussian_check = newLigMOL2Name.replace('.mol2', '.chk')
#
#     # getting atom types and coordinates from the renamed file
#     atoms, coordinate = [], []
#     for line in open(renumbered_pdb):
#         if line.startswith('HETATM') or line.startswith('ATOM'):
#             splitted = line.split()
#             new_line = splitted[2] + ' ' + line[30:54].strip() + '\n'
#             atoms.append(splitted[2][0])
#             coordinate.append(new_line)
#     # writing gaussian input file using the atomtypes and coordinates
#     output = open(gaussian_input, 'w')
#     halogens = {'F': 9, 'Cl': 17, 'Br': 35, 'I': 53, 'At': 85, 'Ts': 117}
#     multiplicity = "1"
#     if any(halogen in ('F', 'Cl', 'Br', 'I', 'At', 'Ts') for halogen in atoms):
#         output.write('%chk=' + gaussian_check + '\n%nprocshared=8\n'
#                                                 '%mem=8GB\n# HF/LANL2DZ Opt=(Redundant) SCF=Tight Geom=PrintInputOrient pop=(mk,ReadRadii) iop(6/33=2,6/41=10,6/42=17)\n'
#                                                 '\n<qmtool> simtype="Geometry optimization" </qmtool>\n'
#                                                 '\n' + f'{wg_total_charge} {multiplicity}\n')
#         for i in coordinate:
#             output.write(i)
#         for atom in atoms:
#             if atom in halogens.keys():
#                 output.write(f'\n{halogens[atom]}  2.0\n')
#     else:
#         output.write('%chk=' + gaussian_check + '\n%nprocshared=8\n'
#                                                 '%mem=8GB\n# HF/6-31G* Opt=(Redundant) SCF=Tight Geom=PrintInputOrient pop=mk iop(6/33=2,6/41=10,6/42=17)\n'
#                                                 '\n<qmtool> simtype="Geometry optimization" </qmtool>\n'
#                                                 '\n' + f'{wg_total_charge} {multiplicity}\n')
#         for i in coordinate:
#             output.write(i)
#         output.write('\n')
#     output.close()
#     # now we run gaussian to prepare the antechamber input file
#     try:
#         run(f'g09 {gaussian_input}', shell=True, check=True)
#     except CalledProcessError as e:
#         print("Gaussian09 failed with error status:", e)
#         run(f'antechamber -i {gaussian_output} -fi gout -o {newLigMOL2Name} -fo mol2 -c resp -eq 2 -s 2 -pf y  -rn UNL;rm esout punch qout QOUT',
#             stdout=DEVNULL)
#
#
# def PrepareProtein():
#     tick = None
#     with open("./receptor/receptor.pdb", 'r') as originalREC:
#         for line in originalREC:
#             if 'POPC' or 'MEMB' or 'POPE' in line:
#                 tick = 32.0
#
#     protein = Molecule('./receptor/receptor.pdb')
#     receptor = systemPrepare(protein, hydrophobic_thickness=tick, ignore_ns_errors=True,
#                              hold_nonpeptidic_bonds=True, titration=True)
#     receptor.write('./receptor/receptor_H.pdb')
#     # we need to remove hydrogens as systemPrepare adds CHARMM-like Hs atomtypes to the system...
#     run('pdb4amber -i ./receptor/receptor_H.pdb -o ./receptor/receptor.pdb -y -a; rm ./receptor/receptor_H.pdb',
#         shell=True)
#     del protein, receptor  # clearing memory as we don't need those anymore
#
#
# def CheckAndRunParameterize() -> None:
#     if not os.path.exists('./receptor/receptor.pdb'):
#         PrepareProtein()
#     if len(ResultsFolders) != 0:
#         with Pool(processes=4) as p:
#             processes = []
#             for foldersLeft in ResultsFolders:
#                 for file in os.listdir(foldersLeft):
#                     if file.endswith('.pdb') and file.startswith(str(foldersLeft).split("/")[-1]) and file.endswith(
#                             'pose1.pdb'):
#                         # trovo in /scratch/ludovico9/galectin/benchmark/post_Docks_amber/ZINC000238950253 = folderLeft
#                         process = p.apply_async(ParameterizeAndCreateComplex, (file, foldersLeft))
#                         processes.append(process)
#             for proc in processes:
#                 proc.get()
#         p.close()
#         p.join()
#
#
# def main():
#     # FindAndMoveLigands()
#     CheckAndRunParameterize()
#     os.chdir(cwd)
#     # torno nel main thread perché multiprocessing è affidabile quanto Salvini
#
#
# if __name__ == '__main__':
#     main()
