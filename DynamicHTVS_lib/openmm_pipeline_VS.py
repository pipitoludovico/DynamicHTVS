import signal
import numpy as np

import openmm as mm
import openmm.app as app
from openmm.app import *
from openmm.unit import *
from openmm import *

import argparse
import subprocess
import logging

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
file_handler = logging.FileHandler('featurizer.log')
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
logging.getLogger().addHandler(file_handler)

try:
    import GPUtil
except ImportError(GPUtil):
    subprocess.Popen('pip install GPUtil', shell=True)


def getpid():
    mypid = os.getpid()
    pidFile = open(".mypid", "w")
    pidFile.write(str(mypid))
    pidFile.close()


def getGPUids(_excluded):
    GPUs = GPUtil.getGPUs()
    gpu_ids = []

    for availableGPU in GPUs:
        gpu_ids.append(availableGPU.id)
    if len(gpu_ids) != 0:
        print("Available GPUS: ", gpu_ids)
    return gpu_ids


def createBatches(b_replicas, b_total_gpu_ids):
    quotient, rest = divmod(b_replicas, len(b_total_gpu_ids))
    b_result = quotient * b_total_gpu_ids + b_total_gpu_ids[:rest]
    b_batches = [b_result[b_i:b_i + len(b_total_gpu_ids)] for b_i in range(0, len(b_result), len(b_total_gpu_ids))]
    return b_batches, b_total_gpu_ids


def runOpenmmEquilibration(eq_coordinates, eq_topology, e_availableIDs, e_userPar, e_userTop, e_protSteps, e_eqSteps,
                           eq_SaveF, protRestraints_, membRestraints_, ligresname_):
    coords, top, charmm, params, amber, gromacs, ligand, proteinCA, membranePOPC, membranePOPE = None, None, None, None, None, None, None, None, None, None
    m_integrator = mm.LangevinMiddleIntegrator(310 * kelvin, 1 / picosecond, 0.002 * picoseconds)
    platform = mm.Platform.getPlatformByName('CUDA')
    # properties = {'DeviceIndex': f'{str(e_availableIDs)}', 'Precision': 'mixed'}
    logging.debug(f"Building Properties...")
    properties = {'DeviceIndex': '0', 'Precision': 'mixed'}
    print('Running thermalization for ', e_protSteps)
    print('Running NPT equilibration for  ', e_eqSteps)
    print("Saving every ", eq_SaveF, " steps.\n")
    if eq_topology.endswith('.prmtop'):
        print('Building from AMBER system...')
        coords = AmberInpcrdFile(eq_coordinates)
        top = AmberPrmtopFile(eq_topology)
        amber = True

    if eq_topology.endswith('.top'):
        print('Building from GROMACS system...')
        coords = GromacsGroFile(eq_coordinates)
        top = GromacsTopFile(eq_topology, periodicBoxVectors=coords.getPeriodicBoxVectors(), includeDir='./oplsaam.ff/')
        gromacs = True

    if any(file_.endswith('.psf') for file_ in os.listdir(os.getcwd())):
        print('Building from CHARMM system...')
        coords = app.PDBFile(eq_coordinates)
        top = app.CharmmPsfFile(eq_topology)
        pos = coords.positions.value_in_unit(nanometers)
        boxLength = np.max(pos, axis=0) - np.min(pos, axis=0)
        x, y, z = boxLength[0], boxLength[1], boxLength[2]
        top.setBox(x * nanometer, y * nanometer, z * nanometer)
        defaultParams = ["/home/scratch/MD_utilities/toppar_c36_jul20/par_all36_prot.prm",
                         "/home/scratch/MD_utilities/toppar_c36_jul20/par_all36_cgenff.prm",
                         "/home/scratch/MD_utilities/toppar_c36_jul20/top_all36_cgenff.rtf",
                         "/home/scratch/MD_utilities/toppar_c36_jul20/par_all36_lipid.prm",
                         "/home/scratch/MD_utilities/toppar_c36_jul20/top_all36_lipid.rtf"]
        if e_userPar is not None and e_userTop is not None:
            user_parameters = [*e_userPar, *e_userTop]
            paramPATHS = defaultParams + user_parameters
        else:
            paramPATHS = defaultParams
        params = app.CharmmParameterSet(*paramPATHS)
        charmm = True

    if charmm:
        system = top.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=0.9 * nanometer,
                                  switchDistance=0.75 * nanometer, constraints=HBonds, rigidWater=True,
                                  hydrogenMass=4 * amu)
    else:
        system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9 * nanometer, switchDistance=0.75 * nanometer,
                                  constraints=HBonds, rigidWater=True, hydrogenMass=4 * amu)
    logging.debug(f"Building System...")

    simulation_eq = Simulation(top.topology, system, m_integrator, platform, properties)

    simulation_eq.context.setPositions(coords.positions)
    logging.debug(f"Context is Set...")

    e_totalSteps = int(e_protSteps + e_eqSteps)

    print("\nCreating the MonteCarlo Barostat...")
    system.addForce(mm.MonteCarloBarostat((1 * bar), (310 * kelvin)))

    def ApplyRestraints(simulation_restraints, restraintIndexesLocal):
        print('\nAdding restraints...')
        restraint = mm.CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
        restraint.addGlobalParameter('k', 10.0 * kilojoules_per_mole / nanometer)
        restraint.addGlobalParameter('t0', 0.0 * picoseconds)
        restraint.addPerParticleParameter('x0')
        restraint.addPerParticleParameter('y0')
        restraint.addPerParticleParameter('z0')
        for atom in restraintIndexesLocal:
            restraint.addParticle(atom, coords.positions[atom].value_in_unit(nanometers))
        system.addForce(restraint)
        simulation_restraints.context.reinitialize(preserveState=True)

    if gromacs:
        if e_userPar:
            print("GROMACS separate ligand file Found. Adding restraints")
            ligPrmtop = []
            for par in e_userPar:
                ligPrmtop.append(par)
        else:
            print("Looking for ligand in the complex topology")

    protResname = ['PRO', 'TRP', 'LYS', 'GLN', 'TYR', 'GLY', 'HIS', 'ARG', 'LEU', 'GLU', 'PHE', 'CYS', 'ILE',
                   'ASN', 'SER', 'THR', 'ASP', 'ALA', 'VAL', 'MET']
    proteinRes = []
    if protRestraints_:
        for protSel in protRestraints_:
            if protSel.split(",")[0] == 'is':
                proteinRes += [p.index for p in top.topology.atoms() if
                               (p.name == protSel.split(",")[1] and p.residue.name in protResname)]
            else:
                proteinRes += [p.index for p in top.topology.atoms() if
                               (p.name != protSel.split(",")[1] and p.residue.name in protResname)]
    membraneRes = []
    if membRestraints_:
        for membSel in membRestraints_:
            print('Adding restraints to:', membSel.split(',')[0], "name", membSel.split(",")[1])
            membraneRes += [p.index for p in top.topology.atoms() if
                            (p.residue.name == membSel.split(",")[0] and p.name == membSel.split(",")[1])]
            print("membraneres", membraneRes)
    ligandRes = []
    if charmm and e_userPar is not None:
        print("CHARMM ligand Found. Adding restraints")
        resNames = []
        if e_userPar and e_userTop:
            for ut in e_userTop:
                with open(ut) as userTopLigandFile:
                    for line in userTopLigandFile.readlines():
                        if 'RESI' in line:
                            resNames.append(line.split()[3])
        ligandRes = [p.index for name in resNames for p in coords.topology.atoms() if name == p.residue.name]
    if amber or gromacs:
        print(f"\nLooking for any {[*ligresname_]} resnames inside the topology...")
        for singleLigResname in ligresname_:
            ligandRes += [p.index for p in top.topology.atoms() if (singleLigResname == p.residue.name)]

    restraintIndexes = list(set([atomIndex for component in (proteinRes, membraneRes, ligandRes) if
                                 component is not None for atomIndex in component]))
    if proteinRes:
        print("Adding restraints to the protein")
    if membraneRes:
        print("Adding restraints to the selected membrane residues/atoms", membraneRes)
    if ligandRes:
        print("Adding restraints to the ligand")
    print("Number of restrained atoms:\n", len(restraintIndexes))
    ApplyRestraints(simulation_eq, restraintIndexes)
    logging.debug(f"Restraints Added...")

    simulation_eq.reporters.append(
        app.StateDataReporter('equilibration.std', 1000, step=True, totalSteps=e_totalSteps, remainingTime=True,
                              potentialEnergy=True, temperature=True))
    simulation_eq.reporters.append(app.DCDReporter('equilibration.dcd', eq_SaveF, enforcePeriodicBox=True))
    logging.debug(f"Reporters Created...")

    print('\nMinimizing local energy...')
    simulation_eq.minimizeEnergy(tolerance=0.1 * kilojoule / mole)
    print('\nWarming up the complex with the restraints...')
    logging.debug(f"Energy Minimized...")

    T = 5
    mdstepsT = 6200
    for t_i in range(1, 63):
        temperature = (t_i * T) * kelvin
        m_integrator.setTemperature(temperature)
        simulation_eq.context.setVelocitiesToTemperature(temperature)
        simulation_eq.step(int(mdstepsT / 62))
    logging.debug(f"Protocol Finished...")

    print('\nEquilibrating with restraints and temperature...')
    easing = int(e_protSteps * 0.2)
    prot = int(e_protSteps - easing)
    simulation_eq.step(prot)
    logging.debug(f"Restrained EQUILIBRATION Finished...")

    print(f'\n{easing} steps easing "k" restraint...:')
    count = 0
    while count < easing:
        for i_T in range(50, -1, -2):
            r_i = (round(i_T * 0.1, 3))
            simulation_eq.context.setParameter('k', r_i * kilocalories_per_mole / nanometer ** 2)
            simulation_eq.step(1000)
            count += 1000
            logging.debug(f"Easing restraints with k= {r_i}")
    logging.debug(f"\nEasing restrained EQUILIBRATION FINISHED")
    print('\nThermalization completed. Equilibrating NPT with k=0 every: ', eq_SaveF, 'steps. Equilibrating for: ',
          e_eqSteps)
    simulation_eq.context.setParameter('k', 0)
    logging.debug(f"\nContext SET")

    simulation_eq.reporters.append(app.CheckpointReporter('equilibration_checkpnt.chk', e_eqSteps))
    # simulation.reporters.append(app.DCDReporter('NPT_K0.dcd', 1000, enforcePeriodicBox=True))
    simulation_eq.step(e_eqSteps)
    logging.debug(f"\nNPT EQUILIBRATION FINISHED")

    # saving the final coordinates
    positions = simulation_eq.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
    app.PDBFile.writeFile(simulation_eq.topology, positions, open("last_frame.pdb", 'w'))
    # saving final state
    final_state = simulation_eq.context.getState(getPositions=True, getVelocities=True)

    with open('equilibration_checkpnt.xml', 'w') as output:
        output.write(XmlSerializer.serialize(final_state))
    logging.debug(f"\nFINAL STATE SAVED")
    simulation_eq.reporters.clear()


def runOpenmmProduction(eq_coordinates, p_topology, p_userPar, p_userTop, p_replicas, p_GPU, m_number_of_steps,
                        m_timeStep, p_saveF):
    os.makedirs(f"run_{p_replicas}", exist_ok=True)
    p_coords, p_top, charmm, params, gromacs = None, None, None, None, None
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'DeviceIndex': f'{str(p_GPU)}', 'Precision': 'mixed'}
    integrator = mm.LangevinIntegrator(310 * kelvin, 1 / picosecond, m_timeStep * picoseconds)
    logging.debug(f"\nPLATFORM INITIALIZED")
    if p_topology.endswith(".top") and eq_coordinates.endswith('.gro'):  # to be tested with gromacs
        gromacs = True
        p_coords = GromacsGroFile(eq_coordinates)
        p_top = GromacsTopFile(p_topology, periodicBoxVectors=p_coords.getPeriodicBoxVectors(),
                               includeDir='./oplsaam.ff/')
    if p_topology.endswith('.prmtop'):
        p_top = AmberPrmtopFile(p_topology)
    if p_topology.endswith('.psf'):
        charmm = True
        p_coords = app.PDBFile(eq_coordinates)
        pos = p_coords.positions.value_in_unit(nanometers)
        boxLength = np.max(pos, axis=0) - np.min(pos, axis=0)
        x, y, z = boxLength[0], boxLength[1], boxLength[2]
        p_top = app.CharmmPsfFile(p_topology)
        p_top.setBox(x * nanometer, y * nanometer, z * nanometer)
        defaultParams = ["/home/scratch/MD_utilities/toppar_c36_jul20/par_all36_prot.prm",
                         "/home/scratch/MD_utilities/toppar_c36_jul20/par_all36_cgenff.prm",
                         "/home/scratch/MD_utilities/toppar_c36_jul20/top_all36_cgenff.rtf",
                         "/home/scratch/MD_utilities/toppar_c36_jul20/par_all36_lipid.prm",
                         "/home/scratch/MD_utilities/toppar_c36_jul20/top_all36_lipid.rtf"]
        if p_userPar is not None and p_userTop is not None:
            user_parameters = [*p_userPar, *p_userTop]
            paramPATHS = defaultParams + user_parameters
        else:
            paramPATHS = defaultParams
        params = app.CharmmParameterSet(*paramPATHS)
    logging.debug(f"\nPRODUCTION PARAMETERS LOADED")

    if gromacs:
        p_top = GromacsTopFile(p_topology, periodicBoxVectors=p_coords.getPeriodicBoxVectors(),
                               includeDir='/usr/local/gromacs/share/gromacs/top')
    if charmm:
        system = p_top.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=0.9 * nanometer,
                                    switchDistance=0.75 * nanometer, constraints=HBonds, rigidWater=True,
                                    hydrogenMass=4 * amu)
    else:
        system = p_top.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9 * nanometer,
                                    switchDistance=0.75 * nanometer, constraints=HBonds, rigidWater=True,
                                    hydrogenMass=4 * amu)
    logging.debug(f"\nSystem BUILT")

    simulation_production = Simulation(p_top.topology, system, integrator, platform, properties)
    simulation_production.loadState("equilibration_checkpnt.xml")
    logging.debug(f"\nPREVIOUS STATE SUCCESSFULLY LOADED")

    simulation_production.reporters.append(
        DCDReporter(f"run_{p_replicas}/Traj_{p_replicas}.dcd", p_saveF, enforcePeriodicBox=True))
    simulation_production.reporters.append(
        StateDataReporter(f'run_{p_replicas}/Traj_{p_replicas}.std', 1000, step=True, totalSteps=m_number_of_steps,
                          remainingTime=True, potentialEnergy=True, temperature=True))
    logging.debug(f"\nPRODUCTION REPORTERS LOADED")

    simulation_production.step(m_number_of_steps)
    logging.debug(f"\nPRODUCTION SIMULATION FINISHED")

    final_state = simulation_production.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    logging.debug(f"\nFINAL STATE CAPTURED.")

    with open(f'run_{p_replicas}/Traj_{p_replicas}.xml', 'w') as output:
        output.write(XmlSerializer.serialize(final_state))
    logging.debug(f"\nFINAL STATE XML SAVED")
    simulation_production.reporters.clear()


def restart(eq_coordinates, p_topology, p_userPar, p_userTop, p_replicas, p_GPU, m_number_of_steps, m_timeStep,
            p_saveF):
    os.makedirs(f"run_{p_replicas}", exist_ok=True)
    p_coords, p_top, charmm, params, gromacs = None, None, None, None, None
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'DeviceIndex': f'{str(p_GPU)}', 'Precision': 'mixed'}
    integrator = mm.LangevinIntegrator(310 * kelvin, 1 / picosecond, m_timeStep * picoseconds)

    if p_topology.endswith(".top") and eq_coordinates.endswith('.gro'):  # to be tested with gromacs
        gromacs = True
        p_coords = GromacsGroFile(eq_coordinates)
        p_top = GromacsTopFile(p_topology, periodicBoxVectors=p_coords.getPeriodicBoxVectors(),
                               includeDir='./oplsaam.ff/')
    if p_topology.endswith('.prmtop'):
        p_top = AmberPrmtopFile(p_topology)
    if p_topology.endswith('.psf'):
        charmm = True
        p_coords = app.PDBFile(eq_coordinates)
        pos = p_coords.positions.value_in_unit(nanometers)
        boxLength = np.max(pos, axis=0) - np.min(pos, axis=0)
        x, y, z = boxLength[0], boxLength[1], boxLength[2]
        p_top = app.CharmmPsfFile(p_topology)
        p_top.setBox(x * nanometer, y * nanometer, z * nanometer)
        defaultParams = ["/home/scratch/MD_utilities/toppar_c36_jul20/par_all36_prot.prm",
                         "/home/scratch/MD_utilities/toppar_c36_jul20/par_all36_cgenff.prm",
                         "/home/scratch/MD_utilities/toppar_c36_jul20/top_all36_cgenff.rtf",
                         "/home/scratch/MD_utilities/toppar_c36_jul20/par_all36_lipid.prm",
                         "/home/scratch/MD_utilities/toppar_c36_jul20/top_all36_lipid.rtf"]
        if p_userPar is not None and p_userTop is not None:
            user_parameters = [*p_userPar, *p_userTop]
            paramPATHS = defaultParams + user_parameters
        else:
            paramPATHS = defaultParams
        params = app.CharmmParameterSet(*paramPATHS)

    if gromacs:
        p_top = GromacsTopFile(p_topology, periodicBoxVectors=p_coords.getPeriodicBoxVectors(),
                               includeDir='/usr/local/gromacs/share/gromacs/top')
    if charmm:
        system = p_top.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=0.9 * nanometer,
                                    switchDistance=0.75 * nanometer, constraints=HBonds, rigidWater=True,
                                    hydrogenMass=4 * amu)
    else:
        system = p_top.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9 * nanometer,
                                    switchDistance=0.75 * nanometer, constraints=HBonds, rigidWater=True,
                                    hydrogenMass=4 * amu)

    simulation_production = Simulation(p_top.topology, system, integrator, platform, properties)
    simulation_production.loadState("equilibration_checkpnt.xml")

    simulation_production.reporters.append(
        DCDReporter(f"run_{p_replicas}/Traj_{p_replicas}.dcd", p_saveF, enforcePeriodicBox=True, append=True))
    simulation_production.reporters.append(
        StateDataReporter(f'run_{p_replicas}/Traj_{p_replicas}.std', 1000, step=True, totalSteps=m_number_of_steps,
                          remainingTime=True, potentialEnergy=True, temperature=True))
    simulation_production.step(m_number_of_steps)

    final_state = simulation_production.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    with open(f'run_{p_replicas}/Traj_{p_replicas}.xml', 'w') as output:
        output.write(XmlSerializer.serialize(final_state))
    simulation_production.reporters.clear()


def GetCoordinatesAndTopology():
    _coordinates, _topology = None, None
    if any(file_.endswith('.psf') for file_ in os.listdir(os.getcwd())):
        _coordinates = [pdb for pdb in os.listdir('./') if pdb.endswith('.pdb')][0]
        _topology = [psf for psf in os.listdir('./') if psf.endswith('.psf')][0]
    if any(file_.endswith('.prmtop') for file_ in os.listdir(os.getcwd())):
        _coordinates = [inpcrd for inpcrd in os.listdir('./') if inpcrd.endswith('.inpcrd')][0]
        _topology = [prmtop for prmtop in os.listdir('./') if prmtop.endswith('.prmtop')][0]
    if any(file_.endswith('.gro') for file_ in os.listdir(os.getcwd())):
        _coordinates = [gro for gro in os.listdir('./') if gro.endswith('.gro')][0]
        _topology = [groTOP for groTOP in os.listdir('./') if groTOP.endswith('.top')][0]
    return _coordinates, _topology


coordinates, topology = GetCoordinatesAndTopology()

print('Coordinate File: ', coordinates)
print('Topology File: ', topology)

ap = argparse.ArgumentParser()
ap.add_argument('-r', '--replicas', type=int, required=False, help=' use -r to set the number of replicas')
ap.add_argument('-gpu', '--gpu', type=str, required=True, help=' use -gpu to set the gpu')
ap.add_argument('-up', '--userpar', type=str, action='append', required=False,
                help=' use -up to set the user parameter')
ap.add_argument('-ut', '--usertop', type=str, action='append', required=False, help=' use -ut to set the user topology')
ap.add_argument('-e', '--excluded', nargs='*', required=False,
                help=' use -e to exclude a list of GPUs from being used by OpenMM: e.g. -e 0 3')
ap.add_argument('-pt', '--protocoltime', type=float, required=False,
                help='set -pt to determine the duration for your thermalization/restraint protocol in ns')
ap.add_argument('-et', '--equilibrationtime', type=float, required=False,
                help='set -et to determine the duration for your equilibration in ns')
ap.add_argument('-prt', '--productiontime', type=float, required=False,
                help='set -prt to determine the duration for your production in ns')
ap.add_argument('-sv', '--savefreq', type=int, required=False,
                help='set -sv to determine the save frequency for your dynamics in ps')
ap.add_argument('-proteinRestraints', '--proteinRestraints', nargs='*', required=False, default=["is,CA", "not,H"],
                help='defines which atoms of the proteins you want to restrain. The synthats is comma-separated. "not,H" or "is,CA". [Default = "is,CA"]')
ap.add_argument('-membraneRestraints', '--membraneRestraints', nargs='*', required=False,
                help='defines the RESNAME and the atom of the membrane you want to restrain. '
                     '\nThe synthats is comma-separated. "POPE,P" or "POPE,P,POPC,P,POPE,N" for multiple selections. [Default = None]')
ap.add_argument('-ligresname', '--ligresname', type=str, action='append', default=["UNL", "UNK"], required=False,
                help=' set -ligresname UNK to identify your small molecule resname [Default = UNL, UNK]')

ap.add_argument('-in', '--integrator', type=float, required=False, help=' use -in to set the timestep in fs')
ap.add_argument("-k", '--kill', required=False, action='store_true', default=False, help="kill the current process")
ap.add_argument('-eq', '--eq', required=False, default=False, help='set -run to run the production')
ap.add_argument('-run', '--run', required=False, default=False, help='set -run to run the production')
ap.add_argument('-restart', '--restart', required=False, default=False,
                help='set -restart to restart the production from the last checkpoint')

args = ap.parse_args()

if args.kill is True:
    import os

    os.system('val=$(<.mypid ) && kill -9 $val')
    os.kill(os.getpid(), signal.SIGKILL)
getpid()

if args.replicas is not None:
    replicas = args.replicas
else:
    replicas = 1

userPar, userTop, excluded, protocolTime, equilibrationTime, productionTime, integratorINfs, saveFreq, protRestraints, membRestraints, gpu = (
    args.userpar, args.usertop, args.excluded, args.protocoltime, args.equilibrationtime, args.productiontime,
    args.integrator, args.savefreq, args.proteinRestraints, args.membraneRestraints, args.gpu)

if excluded:
    if gpu in excluded:
        print("\nWARNING: Conflict between exclusion and gpu. The first available GPU will be used instead.")
        gpuTMP = getGPUids(excluded)
        gpu = int(gpuTMP[0])

if protocolTime is not None:
    protocolSteps = int((protocolTime * 10 ** 6) / 2)
else:
    protocolSteps = int((2 * 10 ** 6) / 2)

if equilibrationTime is not None:
    equilibrationSteps = int((equilibrationTime * 10 ** 6) / 2)
else:
    equilibrationSteps = int(2 * 10 ** 6 / 2)

if integratorINfs is not None:
    integratorINps = integratorINfs / 1000
else:
    integratorINfs = 0.004
    integratorINps = 4 / 1000

if productionTime is not None:
    simulationSteps = int(((productionTime * 10 ** 6) / integratorINfs))
else:
    simulationSteps = int(200 * 10 ** 6 / 4)

if saveFreq is None:
    saveFreq = int(200 / integratorINps)
else:
    saveFreq = int(saveFreq / integratorINps)

print("Settings:")
print("Thermalization Time (in steps): ", protocolSteps)
print("Equilibration Time (in steps): ", equilibrationSteps)
print("Production Time (in steps): ", simulationSteps)
print("Integrator Timestep (in ps): ", integratorINps)
print("Save Frequency (in steps)", saveFreq)
print("User Topology: ", userTop)
print("User Parameters: ", userPar)
print("Protein Restraints: ", protRestraints, "Membrane Restraints: ", membRestraints)

if args.eq:
    runOpenmmEquilibration(coordinates, topology, gpu, userPar, userTop, protocolSteps, equilibrationSteps,
                           saveFreq, protRestraints, membRestraints, args.ligresname)
    print("Equilibration completed.")

if args.run:
    runOpenmmProduction(coordinates, topology, userPar, userTop, 1, gpu, simulationSteps, integratorINps, saveFreq)
    print(f"\nRun finished.")

if args.restart:
    restart(coordinates, topology, userPar, userTop, 1, gpu, simulationSteps, integratorINps, saveFreq)
    print(f"\nRestart finished.")
