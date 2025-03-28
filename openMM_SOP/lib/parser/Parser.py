import argparse
import os.path
from argparse import Namespace
from os import listdir, path, kill, system, getpid
from lib.SOP_utilities.utilityFunctions import getpid, ParseMembraneRestraints
from signal import SIGKILL


def ParseCLIarguments():
    """
    :param: generic types that determine the length of the SOP as well as some options
    :return: dictionary
    """
    ap = argparse.ArgumentParser()

    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("-cyt", "-cytosolic", action='store_true')
    group.add_argument("-memb", "-membrane", action='store_true')

    # user-defined parameters
    group.add_argument("-user", "-userDefined", action='store_true')
    ap.add_argument('-pt', '--protocoltime', type=int, required=False,
                    help='set -pt to determine the duration for your thermalization/restraint protocol in ns')
    ap.add_argument('-et', '--equilibrationtime', type=int, required=False,
                    help='set -et to determine the duration for your equilibration in ns')
    ap.add_argument('-ei', '--equilibrationintegrator', type=float, required=False, default=2.0,
                    help='set -ei to set the integrator ts for your equilibration in fs [Default = 2]')
    ap.add_argument('-prt', '--productiontime', type=int, required=False,
                    help='set -prt to determine the duration for your production in ns')
    ap.add_argument('-pi', '--productionintegrator', type=float, required=False, default=4.0,
                    help='set -pi to set the integrator ts for your production in fs [Default = 4]')

    ap.add_argument('-sv', '--savefreq', type=int, required=False,
                    help='set -sv to determine the save frequency for your dynamics in ps')
    ap.add_argument('-proteinRestraints', '--proteinRestraints', nargs='*', required=False, default=["is,CA"],
                    help='defines which atoms of the proteins you want to restrain. The synthats is comma-separated. "not,H" or "is,CA". [Default = "is,CA"]')
    ap.add_argument("-membraneRestraints", "--membraneRestraints", nargs="+", action=ParseMembraneRestraints,
                    help='defines the RESNAME and the atom of the membrane you want to restrain.'
                         '\nThe synthats is comma-separated. "POPE,P" or "POPE,P,POPC,P,POPE,N" for multiple selections. [Default = None]')
    ap.add_argument('-ligresname', '--ligresname', type=str, action='append', default=None, required=False,
                    help=' set -ligresname UNK to identify your small molecule resname [Default = UNL, UNK, LIG]')
    ap.add_argument("-FEP", action='store_true', required=False,
                    help="Attempt a direct calculation of the alchemical free energy for a ligand.")
    ap.add_argument("-eqBarostat", '--eb', required=False, action='store_false',
                    help="applies Montecarlo barostat to the equilibration protocol")
    ap.add_argument("-prBarostat", '--pb', required=False, action='store_true',
                    help="applies Montecarlo barostat to the production protocol")

    ap.add_argument('-r', '--replicas', type=int, required=False, default=1,
                    help=' use -r to set the number of replicas')
    ap.add_argument('-ex', '--excluded', nargs='*', required=False,
                    help=' use -e to exclude a list of GPUs from being used by OpenMM: e.g. -e 0 3')
    ap.add_argument("-k", '--kill', required=False, action='store_true', default=False, help="kill the current process")
    ap.add_argument('-eq', '--eq', required=False, action='store_true', help='set -eq to run the equilibrating SOP')
    ap.add_argument('-rEq', '--recoverEquilibration', required=False, action='store_true',
                    help='set -rEq to recover the equilibration from the last step')
    ap.add_argument('-exEq', '--extendEquilibration', required=False, type=int,
                    help='set -exEq to extend the equilibration')
    ap.add_argument('-run', '--run', required=False, action='store_true', help='set -run to run the production SOP')
    ap.add_argument('-rPr', '--recoverProduction', required=False, action='store_true',
                    help='set -rPr to recover the production from the last step')
    ap.add_argument('-exPr', '--extendProduction', required=False, type=int, help='set -exPr to extend the production')
    ap.add_argument('-GPU', '--GPU', required=False, type=str, help='set -GPU to specify a single GPU id to use')

    args: Namespace = ap.parse_args()

    if args.kill is True:
        system('val=$(<.mypid ) && kill -9 $val')
        kill(os.getpid(), SIGKILL)

    getpid()
    protocolTime: int = args.protocoltime if args.protocoltime else 8 if args.cyt else 120 if args.memb else args.protocoltime
    equilibrationTime: int = args.equilibrationtime if args.equilibrationtime else 2 if args.cyt else 30 if args.memb else args.equilibrationtime
    productionTime: int = args.productiontime if args.productiontime else 500 if args.cyt else 1000 if args.memb else args.productiontime
    saveFreq: int = 50000
    proteinRestraints: list = ["is,CA"]
    membraneRestraints: list = [None][0] if args.cyt else ['is,P', 'is,C11'] if args.memb else args.membraneRestraints
    ligResname: list = ["UNL", "UNK", 'LIG'] if args.ligresname is None else args.ligresname
    equilibrationBarostat: bool = args.eb
    productionBarostat: bool = args.pb
    equilibrationIntegrator: int = args.equilibrationintegrator if args.equilibrationintegrator else 2
    productionIntegrator: int = args.productionintegrator if args.productionintegrator else 4

    arguments = {"replicas": args.replicas, 'excludedGPUs': args.excluded, 'protocolTime': protocolTime,
                 'eqBarostat': equilibrationBarostat, 'prBarostat': productionBarostat,
                 'equilibrationTime': equilibrationTime, 'productionTime': productionTime,
                 'eqIntegrator': equilibrationIntegrator, 'prIntegrator': productionIntegrator,
                 'saveFreq': saveFreq, 'proteinRestraints': proteinRestraints,
                 'membraneRestraints': membraneRestraints, 'ligandResname': ligResname, 'eq': args.eq,
                 'run': args.run, 'recoverEq': args.recoverEquilibration, 'recoverPr': args.recoverProduction,
                 'extendEq': args.extendEquilibration, 'extendPr': args.extendProduction,
                 'GPU': args.GPU, 'FEP': args.FEP}
    return arguments


def ParseFiles():
    coordinates, topology = None, None
    userTopologies, userParameters = [], []
    topologyExtensions = ('psf', 'prmtop', 'top')
    coordinatesExtensions = ('pdb', 'gro', 'inpcrd')
    for file in listdir('.'):
        if file.endswith(topologyExtensions):
            topology = path.abspath(file)
        if file.endswith(coordinatesExtensions):
            coordinates = path.abspath(file)
    if path.exists('./extraTopPar'):
        extraTopExtensions = ('rtf', 'str', 'top')
        extraParExtensions = ('par', 'prm', 'top')
        for extraFile in listdir('./extraTopPar'):
            if extraFile.endswith(extraTopExtensions):
                userTopologies.append(os.path.abspath('./extraTopPar/' + extraFile))
            if extraFile.endswith(extraParExtensions):
                userParameters.append(os.path.abspath('./extraTopPar/' + extraFile))

    fileDict = {'coordinates': coordinates, 'topology': topology,
                'extraTopologies': userTopologies, 'extraParameters': userParameters}
    print(fileDict)
    return fileDict
