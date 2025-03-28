from multiprocessing.pool import ThreadPool

from lib.parser.Parser import ParseCLIarguments, ParseFiles
from lib.builder.systemBuilder import BuildSystem
from lib.SOP_utilities.utilityFunctions import DescribeSettings, getGPUids, createBatches
from lib.runner.equilibration import Equilibrate
from lib.runner.production import RunProduction


def main():
    parameters = ParseCLIarguments()
    files = ParseFiles()
    DescribeSettings(arguments=parameters)
    system, top, coords, FF = BuildSystem(coordinates=files['coordinates'], topology=files['topology'],
                                          userTop=files['extraTopologies'], userPar=files['extraParameters'])
    availableIDs = getGPUids(parameters['excludedGPUs'])

    if parameters['eq'] or parameters['recoverEq'] or parameters['extendEq']:
        Equilibrate(system, top, coords, files, FF, parameters, availableIDs)

    if parameters['run'] or parameters['recoverPr'] or parameters['extendPr']:
        GPUbatches, idList = createBatches(parameters['replicas'], availableIDs)
        i = 0
        with ThreadPool(processes=len(idList)) as pool:
            results = []
            for GPUbatch in GPUbatches:
                for GPU in GPUbatch:
                    results.append(pool.apply_async(RunProduction, args=(GPU, i, system, top, parameters, FF, files,coords)))
                    i += 1
            for result in results:
                result.get()
        pool.join()
        pool.close()
        pool.terminate()
        print(f"All batches finished.")


if __name__ == '__main__':
    main()
