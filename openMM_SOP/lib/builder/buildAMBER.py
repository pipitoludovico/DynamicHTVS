from openmm.unit import *
from openmm.app import AmberInpcrdFile, AmberPrmtopFile
from openmm.app import PME, HBonds


def amberBuilder(eq_coordinates, eq_topology):
    print('Building from AMBER system...')
    coords = AmberInpcrdFile(eq_coordinates)
    top = AmberPrmtopFile(eq_topology)
    system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9 * nanometer, switchDistance=0.75 * nanometer,
                              constraints=HBonds, rigidWater=True, hydrogenMass=4 * amu)
    return system, top, coords
