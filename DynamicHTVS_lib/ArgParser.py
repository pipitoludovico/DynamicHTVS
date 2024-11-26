import argparse
import shutil
from os import remove
import GPUtil

from .Utility import *


class ArgParser:
    def __init__(self):
        GetRightSettings()
        self.argumentParser()

    @staticmethod
    def argumentParser():
        ap = argparse.ArgumentParser()
        ap.add_argument("-d", '--dock', action='store_true', required=False,
                        help="set -d True to rerun the docking [Default=False] ")
        ap.add_argument("-r", '--rank', action='store_true', required=False,
                        help="set -r True to rerun the ranking [Default=False] ")
        ap.add_argument("-p", '--parameterize', action='store_true', required=False,
                        help="set -p True to rerun the ligand parametrization [Default=False]")
        ap.add_argument("-b", '--build', action='store_true', required=False,
                        help="set -b True to rerun the system builder and MDs [Default=False] ")
        ap.add_argument("-ligands", '--ligands', type=dirOrfile, required=False, default=('smi', 'smi'),
                        help="set -pdb to mark a folder with pdbs or a single pdb [Default=False] ")
        ap.add_argument("-launch", "--launchOpenmm", action='store_true', required=False,
                        help="add -launch to run OpenMM with the prepared files. [Default=False]")
        ap.add_argument("-c", '--calculatescore', action='store_true', required=False,
                        help="set -c True to rerun the Dynamic Score calculation [Default=False] ")
        ap.add_argument("-amber", '--amber', action='store_true', required=False,
                        help="set -amber to use AMBER instead of CHARMM [Default=False] ")
        ap.add_argument("-batch", "-batch", required=False, type=int, default=None, choices=range(1, len(GPUtil.getGPUs()) + 1),
                        help="add -batch to specify the number of parallel simulations you want to run on your machine [Default == number of available GPUS]")
        ap.add_argument("-consider", "-consider", required=False, type=int, default=20,
                        help="add -consider to specify the number compounds you want to consider for MD from your docking results [Default == 20]")
        ap.add_argument("-boxsize", "-bs", required=False, type=int, default=12,
                        help="add -boxsize -bs to set the grid box size [Default == 12]")
        ap.add_argument("-poses", "-poses", required=False, type=int, default=1, help="add -poses to set how many docking poses you want to generate")
        ap.add_argument('-e', '--exclude', nargs='*', required=False, help=' use -e to exclude a list of GPUs from being used by the Dynamic Screener: e.g. -e 0 3')
        ap.add_argument("-sel", '--selection', required=False, type=str, nargs="+",
                        help="a selection -sel between quotes is required for the docking using VMD's selection language.")
        ap.add_argument("-prt", "-prt", required=False, type=int, default=10,
                        help='Set for how long you want to run the post-docking dynamics in ns [Default == 10]')
        ap.add_argument("-purge", "--purge", action='store_true',
                        help='add -purge to remove all the results and start from scratch. WARNING: THIS WILL DELETE ALL THE CONTENT OF THE FOLDER EXEPT FOR THE .smi FILE and equilibration folder!!!')
        ap.add_argument('-membraneRestraints', '--membraneRestraints', nargs='*', required=False,
                        help='defines the RESNAME and the atom of the membrane you want to restrain. The synthats is comma-separated. "POPE,P" or "POPE,P,POPC,P,POPE,N" for multiple selections. [Default = None]')

        args = ap.parse_args()
        if args.purge:
            confirmation = input(
                "Do you really want to purge this folder? It will keep only the .smi file and the 'equilibration' folder in this working directory! [Y/N]")
            if confirmation.upper().startswith('Y'):
                for element in listdir("./"):
                    if not path.isdir(element):
                        if not element.endswith('.smi'):
                            remove(element)
                    if path.isdir(element):
                        if element != 'equilibration':
                            shutil.rmtree(element, ignore_errors=True)
            exit()

        dock, rank, parametrize, build, calculate, selection, amber, launch, batch, exclude, prt, ligands, poses = \
            args.dock, args.rank, args.parameterize, args.build, args.calculatescore, args.selection, args.amber, args.launchOpenmm, args.batch, args.exclude, args.prt, args.ligands, args.poses
        if selection:
            joinedSel = ' '.join(args.selection)
        else:
            joinedSel = None

        return dock, rank, parametrize, build, calculate, joinedSel, amber, launch, batch, exclude, args.consider, prt, args.membraneRestraints, args.boxsize, ligands, poses
