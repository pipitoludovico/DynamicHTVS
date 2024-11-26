#!/usr/bin/env python3
import os
from time import perf_counter

from DynamicHTVS_lib import DatabaseDocker
from DynamicHTVS_lib import DockingRanker
from DynamicHTVS_lib import ArgParser
from DynamicHTVS_lib import Utility

import importlib.resources

package_dir = str(importlib.resources.files('DynamicHTVS_lib'))
OPENMM_SCRIPT_PATH = os.path.join(package_dir, 'openmm_pipeline_VS.py')
Dynamic = os.path.join(package_dir, 'DynamicScorer.py')

ap = ArgParser.ArgParser()
dock, rank, parametrize, build, calculate, selection, amber, launch, batch, exclude, consider, prt, membraneList, boxsize, ligands, poses = ap.argumentParser()


def main():
    os.makedirs('logs', exist_ok=True)
    ROOT = os.getcwd()
    if ligands:
        ligandType, fullLigandPath = ligands[0], ligands[1]
    else:
        ligandType, fullLigandPath = None, None

    if dock:
        print("*" * 200)
        # Dock your database
        if not selection:
            print("Please add a selection using the -sel argument")
            exit()
        dbDocker = DatabaseDocker.DatabaseDocker(amber, boxsize, ligandType, fullLigandPath, poses)
        dbDocker.DockMols(selection)
        print("\n\nDocking completed.")
        print("*" * 200)
    os.chdir(ROOT)

    if rank:
        print("*" * 200)
        print("Ranking results...")
        # Rank your results
        receptorPath = DatabaseDocker.DatabaseDocker(amber, boxsize).GetReceptorPath()
        ranker = DockingRanker.DockingRanker()
        ranker.ReadVinaScores()
        ranker.GetContacts(receptorPath)
        ranker.CreateSummaryChart()
        ranker.SortSummary(consider)
        print("\nRanking completed.\n")
        print("*" * 200)
    os.chdir(ROOT)

    #  receptor
    if amber:
        import DynamicHTVS_lib.ReceptorTools.PrepareReceptor_AMBER
        DynamicHTVS_lib.ReceptorTools.PrepareReceptor_AMBER.PrepareProtein()
    else:
        import DynamicHTVS_lib.ReceptorTools.PrepareReceptor_CHARMM
        DynamicHTVS_lib.ReceptorTools.PrepareReceptor_CHARMM.PrepareProtein()

    #  ligands
    if parametrize:
        Result_Folders = Utility.FindAndMoveLigands(amber, consider)
        print("*" * 200)
        if amber:
            print("Building AMBER Ligand Parameters\n\n")
            from DynamicHTVS_lib.LigandTools.Parameterizer_AMBER import RunParameterize
            RunParameterize(Result_Folders)
        else:
            print("Building CHARMM Ligand Parameters\n\n")
            from DynamicHTVS_lib.LigandTools.Parameterizer_CHARMM import RunParameterize_CHARMM
            RunParameterize_CHARMM(Result_Folders)
        print("\nParameterization completed.")

    if build:
        Result_Folders = Utility.GetReultFolders(amber)
        print("@" * 200)
        print("Building systems...")
        if len(Result_Folders) == 0:
            exit(
                "No post_Docks_amber folder found. Make sure you have run the docking first or create this folder manually.")
        if amber:
            from DynamicHTVS_lib.ComplexTools.AMBER_COMPLEX_BUILDER import BuildAMBERsystems
            BuildAMBERsystems(Result_Folders)
        else:
            from DynamicHTVS_lib.ComplexTools.CHARMM_COMPLEX_BUILDER import BuildCHARMMsystems
            BuildCHARMMsystems(Result_Folders)

    if launch:
        Result_Folders = Utility.GetReultFolders(amber)
        from DynamicHTVS_lib.Runners.BUILD_AND_RUN import RunOpenMMbatches
        RunOpenMMbatches(Result_Folders, OPENMM_SCRIPT_PATH, amber, batch, exclude, prt, membraneList)
        print("\nMD completed")

    os.chdir(ROOT)
    # Calculate the Dynamic Score
    if calculate:
        print("")
        print("#" * 200)
        print("Calculating Scores...")
        os.system(f'python {Dynamic} {amber}')
        print("Scoring completed")


if __name__ == "__main__":
    start = perf_counter()
    main()
    end = perf_counter()
    final = end - start
    print("The process took: ", final)
