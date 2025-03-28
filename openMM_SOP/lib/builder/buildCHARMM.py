import numpy as np
from openmm.unit import *
from openmm.app import *
import warnings

warnings.filterwarnings('ignore')


def charmmBuilder(eq_coordinates, eq_topology, e_userPar=None, e_userTop=None):
    print('Building from CHARMM system...')
    coords = PDBFile(eq_coordinates)
    top = CharmmPsfFile(eq_topology)
    pos = coords.positions.value_in_unit(nanometers)
    boxLength = np.max(pos, axis=0) - np.min(pos, axis=0)
    x, y, z = boxLength[0], boxLength[1], boxLength[2]
    top.setBox(x * nanometer, y * nanometer, z * nanometer)
    defaultParams = [
        "/home/scratch/MD_utilities/toppar_c36_jul20/top_all36_prot.rtf",
        "/home/scratch/ludovico/par_all36_prot_NOcomm.prm",
        '/home/scratch/MD_utilities/toppar_c36_jul20/top_all36_na.rtf',
        '/home/scratch/ludovico/par_all36_na.prm',
        "/home/scratch/MD_utilities/toppar_c36_jul20/top_all36_lipid.rtf",
        "/home/scratch/MD_utilities/toppar_c36_jul20/par_all36_lipid.prm",

        "/home/scratch/MD_utilities/toppar_c36_jul20/stream/lipid/toppar_all36_lipid_inositol.str",
        "/home/scratch/MD_utilities/toppar_c36_jul20/stream/carb/toppar_all36_carb_glycolipid.str",
        "/home/scratch/MD_utilities/toppar_c36_jul20/stream/lipid/toppar_all36_lipid_cholesterol.str",
        "/home/scratch/MD_utilities/toppar_c36_jul20/stream/lipid/toppar_all36_lipid_sphingo.str",

        "/home/scratch/MD_utilities/toppar_c36_jul20/par_all36_carb.prm",

        "/home/scratch/MD_utilities/toppar_c36_jul20/par_all36_cgenff.prm",
        "/home/scratch/MD_utilities/toppar_c36_jul20/top_all36_cgenff.rtf",

        "/home/scratch/MD_utilities/toppar_c36_jul20/toppar_water_ions.str",
    ]

    if e_userPar is not None and e_userTop is not None:
        user_parameters = [*e_userPar, *e_userTop]
        paramPATHS = defaultParams + user_parameters
    else:
        paramPATHS = defaultParams
    params = CharmmParameterSet(*paramPATHS)
    system = top.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=0.9 * nanometer,
                              switchDistance=0.75 * nanometer, constraints=HBonds, rigidWater=True,
                              hydrogenMass=4 * amu)
    return system, top, coords
