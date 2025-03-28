from openmm.unit import *
from openmm.app import *


def gromacsBuilder(eq_coordinates, eq_topology):
    print('Building from GROMACS system...')
    coords = GromacsGroFile(eq_coordinates)
    top = GromacsTopFile(eq_topology, periodicBoxVectors=coords.getPeriodicBoxVectors(), includeDir='./oplsaam.ff/')
    system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9 * nanometer, switchDistance=0.75 * nanometer,
                              constraints=HBonds, rigidWater=True, hydrogenMass=4 * amu)
    return system, top, coords
