from openmm import *
import openmm as mm
import numpy as np
from pymbar import MBAR, timeseries

from ..builder.mdSettings import ApplyRestraints
nsteps = 10000  # number of steps per sample
niterations = 5  # number of samples to collect per alchemical state
lambdas = np.linspace(0.0, 1.0, 10)  # alchemical lambda schedule
nstates = len(lambdas)
u_kln = np.zeros([nstates, nstates, niterations], np.float64)


def SetAlchemicalForces(system, ligandName, top, FF, files, coords):
    forces = {force.__class__.__name__: force for force in system.getForces()}
    nbforce = forces['NonbondedForce']

    # Add a CustomNonbondedForce to handle only alchemically-modified interactions
    ligandRes = []
    chemicalEnv = set([p.index for p in top.topology.atoms()])

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

    print(ligandRes)
    ApplyRestraints(system, coords, ligandRes, 10)

    alchemical_particles = set([atomIndex for atomIndex in ligandRes])
    chemical_particles = chemicalEnv - alchemical_particles
    print(f"Alchemical particles idxs: {alchemical_particles}")

    energy_function = '((1-is_alchemical) + lambda*is_alchemical)*4*epsilon*x*(x-1.0);'
    energy_function += 'x = (sigma/reff_sterics)^6; reff_sterics = sigma*(0.5*(1.0-lambda)*is_alchemical + (r/sigma)^6)^(1/6);'
    energy_function += 'sigma = sqrt(sigma1*sigma2); epsilon = sqrt(epsilon1*epsilon2); is_alchemical = 1 - (solv1*solv2);'

    custom_force = openmm.CustomNonbondedForce(energy_function)
    custom_force.addGlobalParameter('lambda', 0.0)
    custom_force.addPerParticleParameter('sigma')
    custom_force.addPerParticleParameter('epsilon')
    custom_force.addPerParticleParameter('solv')

    for index in range(system.getNumParticles()):
        [charge, sigma, epsilon] = nbforce.getParticleParameters(index)
        if index in alchemical_particles:
            # set charge epsilon of Nonbondedforce to 0
            nbforce.setParticleParameters(index, charge * 0.0, sigma, epsilon * 0.0)
            custom_force.addParticle([sigma, epsilon, 0.0])
        elif index in chemical_particles:
            # set charge epsilon of Nonbondedforce to 0
            custom_force.addParticle([sigma, epsilon, 1.0])
            nbforce.setParticleParameters(index, charge, sigma, epsilon * 0.0)

    for i in range(nbforce.getNumExceptions()):
        (p1, p2, q, sig, eps) = nbforce.getExceptionParameters(i)
        # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED FORCE
        custom_force.addExclusion(p1, p2)
    system.addForce(custom_force)

    for force in system.getForces():
        if isinstance(force, mm.CustomNonbondedForce):
            force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            force.setCutoffDistance(10 * unit.angstroms)
            # force.setUseLongRangeCorrection(True)
        elif isinstance(force, mm.NonbondedForce):
            force.setNonbondedMethod(mm.NonbondedForce.PME)
            force.setCutoffDistance(10 * unit.angstroms)
            force.setUseDispersionCorrection(True)
    return system


def RunFEP(alchemic_simulation, integrator):
    simulation = alchemic_simulation
    kT = unit.AVOGADRO_CONSTANT_NA * unit.BOLTZMANN_CONSTANT_kB * integrator.getTemperature()
    simulation.context.reinitialize()
    print("simulation reinitialized")
    for k in range(nstates):
        for iteration in range(niterations):
            print('state %5d iteration %5d / %5d' % (k, iteration, niterations))
            # Set alchemical state
            simulation.context.setParameter('lambda', lambdas[k])
            # Run some dynamics
            simulation.step(nsteps)
            # Compute energies at all alchemical states
            for state_ in range(nstates):
                simulation.context.setParameter('lambda', lambdas[state_])
                state = simulation.context.getState(getEnergy=True, getPositions=True)
                energy = state.getPotentialEnergy()
                volume = state.getPeriodicBoxVolume()
                pressure = 1 * unit.bar
                u_kln[k, state_, iteration] = (energy + pressure * volume * unit.AVOGADRO_CONSTANT_NA) / kT


def CalculateEnergy():
    N_k = np.zeros([nstates], np.int32)
    for k in range(nstates):
        [_, g, _] = timeseries.detect_equilibration(u_kln[k, k, :])
        indices = timeseries.subsample_correlated_data(u_kln[k, k, :], g=g)
        N_k[k] = len(indices)
        u_kln[k, :, 0:N_k[k]] = u_kln[k, :, indices].T
    mbar = MBAR(u_kln, N_k)
    results = mbar.compute_free_energy_differences()
    with open("FEP.txt", 'w') as FEPout:
        FEPout.write(f'dG: {results["Delta_f"][nstates - 1][0]}, ddG: {results["dDelta_f"][nstates - 1][0]}')
    print(f'dG: {results["Delta_f"][nstates - 1][0]}, ddG: {results["dDelta_f"][nstates - 1][0]}')
