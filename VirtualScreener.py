#!/usr/bin/env python3
import os
from time import perf_counter
import importlib.resources

from DynamicHTVS_lib.CLI_Parser.Commands import *
from DynamicHTVS_lib.Docking import DatabaseDocker, DockingRanker
from DynamicHTVS_lib.Utilities import Utility
from DynamicHTVS_lib.Dynamics.Plotter import DataPlotter

package_dir = str(importlib.resources.files('DynamicHTVS_lib'))
OPENMM_SCRIPT_PATH = os.path.join(package_dir, 'Dynamics/openmm_pipeline_VS.py')
Dynamic = os.path.join(package_dir, 'Dynamics/DynamicScorer.py')


def main():
    print("\n")
    print("*" * 50)
    print("*" * 50)
    print("Weighted Ensemble Dynamic Virtual Screening")
    print("*" * 50)
    print("*" * 50)

    os.makedirs('logs', exist_ok=True)
    ROOT = os.getcwd()

    if dock:
        print("*" * 50)
        # Dock your database
        if not selection_:
            print("Please add a selection using the -sel argument")
            exit()
        dbDocker = DatabaseDocker.DatabaseDocker(amber, boxsize, ligands[0], ligands[1],
                                                 poses)  # ligand[0] = type ligand[1] path
        dbDocker.DockMols(selection_)
        print("\n\nDocking completed.")
        print("*" * 50)
    os.chdir(ROOT)

    if rank:
        print("*" * 50)
        print("Ranking results...")
        # Rank your results
        receptorPath = DatabaseDocker.DatabaseDocker(amber, boxsize).GetReceptorPath()
        ranker = DockingRanker.DockingRanker(restrict)
        ranker.ReadVinaScores()
        ranker.GetContacts(receptorPath)
        ranker.CreateSummaryChart()
        ranker.SortSummary(consider)
        print("\nRanking completed.\n")
        print("*" * 50)
    os.chdir(ROOT)

    #  ligands
    if parameterize:
        Result_Folders = Utility.FindAndMoveLigands(amber, consider)
        print("*" * 50)
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
        #  receptor
        if amber:
            import DynamicHTVS_lib.ReceptorTools.PrepareReceptor_AMBER
            DynamicHTVS_lib.ReceptorTools.PrepareReceptor_AMBER.PrepareProtein()
        else:
            import DynamicHTVS_lib.ReceptorTools.PrepareReceptor_CHARMM
            DynamicHTVS_lib.ReceptorTools.PrepareReceptor_CHARMM.PrepareProtein()
        print("@" * 50)
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

    if launchOpenmm:
        Result_Folders = Utility.GetReultFolders(amber)
        from DynamicHTVS_lib.Runners.BUILD_AND_RUN import RunOpenMMbatches
        RunOpenMMbatches(Result_Folders, OPENMM_SCRIPT_PATH, amber, batch, exclude, productionTime, membraneRestraints)
        print("\nMD completed")

    os.chdir(ROOT)
    # Calculate the Dynamic Score
    if calculateScore:
        print("")
        print("#" * 50)
        print("Calculating Scores...")
        os.system(f'python {Dynamic} {amber}')
        print("Scoring completed")
        print("Plotting...")
        dp = DataPlotter(batch=consider)
        dp.GetResults()
        dp.ParallelRMSDs()
        dp.PlotAll()


if __name__ == "__main__":
    start = perf_counter()
    main()
    end = perf_counter()
    final = end - start
    print("The process took: ", final)
