import buildAMBER, buildGROMACS, buildCHARMM


# import lib.builder.buildAMBER
# import lib.builder.buildGROMACS
# import lib.builder.buildCHARMM


def BuildSystem(coordinates, topology, **kwargs):
    """ Generates the system and returns the FF format for building the simulation"""
    system, top, coords, FF = None, None, None, None
    userTop, userPars = kwargs['userTop'], kwargs['userPar']
    if coordinates is None or topology is None:
        exit('Please put your coordinates and topology files in the folder.')
    if topology.endswith('.prmtop'):
        system, top, coords = buildAMBER.amberBuilder(eq_coordinates=coordinates, eq_topology=topology)
        FF = 'AMBER'
    if topology.endswith('.top'):
        system, top, coords = buildGROMACS.gromacsBuilder(eq_coordinates=coordinates, eq_topology=topology)
        FF = 'GROMACS'
    if topology.endswith('.psf'):
        system, top, coords = buildCHARMM.charmmBuilder(eq_coordinates=coordinates, eq_topology=topology,
                                                        e_userPar=userPars,
                                                        e_userTop=userTop)
        FF = 'CHARMM'
    return system, top, coords, FF
