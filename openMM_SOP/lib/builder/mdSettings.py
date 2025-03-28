import openmm as mm
from openmm.unit import *


def GetMDsettings(GPU: str = 0, friction: float = 0.002):
    m_integrator = mm.LangevinMiddleIntegrator(310 * kelvin, 1 / picosecond, friction * picoseconds)
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'DeviceIndex': f"{GPU}", 'Precision': 'mixed'}
    return m_integrator, platform, properties


def SetBarostat(system):
    print("\nCreating the MonteCarlo Barostat...")
    system.addForce(mm.MonteCarloBarostat((1 * bar), (310 * kelvin)))


def RemoveRestraints(system):
    """Remove all CustomExternalForce instances from the system."""
    forces_to_remove = []
    for i, force in enumerate(system.getForces()):
        if isinstance(force, mm.CustomExternalForce):
            forces_to_remove.append(i)
    # Remove in reverse to avoid index shifting
    for i in reversed(forces_to_remove):
        system.removeForce(i)
    return system


def ApplyRestraints(system, coords, restraintIndexesLocal, kcal=5.0):
    print('\nAdding restraints...')
    restraint = mm.CustomExternalForce('k*periodicdistance(x, y, z, x0, y0, z0)^2')
    # unit.kilocalories_per_mole/unit.nanometers**2
    restraint.addGlobalParameter('k', kcal * kilocalories_per_mole / nanometers**2)
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')
    for atom in restraintIndexesLocal:
        restraint.addParticle(atom, coords.positions[atom].value_in_unit(nanometers))
    system.addForce(restraint)


def GetRestraintIndex(protRestraints, membRestraints, ligandName, top, FF, files, coords):
    protResname = ['PRO', 'TRP', 'LYS', 'GLN', 'TYR', 'GLY', 'HIS', 'ARG', 'LEU', 'GLU', 'PHE', 'CYS', 'ILE',
                   'ASN', 'SER', 'THR', 'ASP', 'ALA', 'VAL', 'MET', 'HSP', 'HSD', 'HSE', 'ALAD', 'ALY', "DORN",
                   "NLE", "AIBD", "AIB", "ASPM", "LEM", "LYM", "NMGLYD", "MGY", "HYP", "HIC", "HICP", "NAPH", "XCP",
                   "CYM", "BPA", "CHA", "SEM", "WF1", "CYSP"]
    membraneResnames = ["POPC", "POPE", "POPI", "CHL1", "SOPE", "POPS", "SSM"]
    proteinRes = []

    if protRestraints:
        for protSel in protRestraints:
            if protSel.split(",")[0] == 'is':
                proteinRes += [p.index for p in top.topology.atoms() if
                               (p.name == protSel.split(",")[1] and p.residue.name in protResname)]
            else:
                proteinRes += [p.index for p in top.topology.atoms() if
                               (p.name != protSel.split(",")[1] and p.residue.name in protResname)]
    if len(proteinRes) > 0:
        print("Number of restrained protein atoms: ", len(proteinRes))

    membraneRes = []
    if membRestraints:
        for membSel in membRestraints:
            print('Adding membrane restraints to:', "atom name: ", membSel.split(",")[1])
            membraneRes += [p.index for p in top.topology.atoms() if (p.residue.name in membraneResnames and p.name == membSel.split(",")[1]) or (p.residue.chain.id == "M" and p.name == membSel.split(",")[1])]
    if len(membraneRes) > 0:
        print("Number of restrained membrane atoms: ", len(proteinRes))
        with open("restrained_membrane_atoms.txt", 'w') as membList:
            for x in membraneRes:
                membList.write(f"{x}\n")

    ligandRes = []
    if FF == "CHARMM":
        resNames = []
        extraTopFile = files.get('extraTopologies')
        for file in extraTopFile:
            with open(file, 'r') as userTopLigandFile:
                for line in userTopLigandFile.readlines():
                    if 'RESI' in line:
                        lig = str(line.split()[1])
                        resNames.append(lig)
                        print("\nCHARMM ligand Found. Adding restraints to resname: ", lig, "\n")
        ligandRes = [p.index for name in resNames for p in coords.topology.atoms() if name == p.residue.name]

    print(f"\nLooking for any {[*ligandName]} resnames inside the topology...")
    size_ligandIdx = 0
    for singleLigResname in ligandName:
        ligandRes += [p.index for p in top.topology.atoms() if (singleLigResname == p.residue.name)]
        if len(ligandRes) > size_ligandIdx:
            size_ligandIdx += len(ligandRes)
            print(f"Found a ligand with resname {singleLigResname}. Adding indexes for the harmonic restraining.")
    restraintIndexes = list(set([atomIndex for component in (proteinRes, membraneRes, ligandRes) if component is not None for atomIndex in component]))

    print("Number of restrained ligand atoms: ", len(ligandRes))
    return restraintIndexes
