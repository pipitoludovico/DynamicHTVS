import os.path

from ..builder.mdSettings import GetMDsettings, SetBarostat, GetRestraintIndex, ApplyRestraints

from openmm.app import *
from openmm import XmlSerializer
from openmm.unit import *


def Equilibrate(system, top, coords, files, FF, parameters, availableIDs) -> None:
    protocolSteps = int((parameters['protocolTime'] * 10 ** 6) / parameters['eqIntegrator'])
    equilibrationSteps = int((parameters['equilibrationTime'] * 10 ** 6) / parameters['eqIntegrator'])
    e_totalSteps = int(protocolSteps + equilibrationSteps)
    allIDs = list(map(str, availableIDs))
    allIDs = ",".join(allIDs)

    def InstantiateSystem():
        friction: float = (parameters['eqIntegrator']) / 1000
        m_integrator, platform, properties = GetMDsettings(parameters.get('GPU', allIDs), friction)
        if parameters['eqBarostat']:
            SetBarostat(system)
            print("\nAdded MonteCarlo Barostat")
        restraintIdxs = GetRestraintIndex(protRestraints=parameters['proteinRestraints'],
                                          membRestraints=parameters['membraneRestraints'],
                                          ligandName=parameters['ligandResname'], top=top, FF=FF, files=files,
                                          coords=coords)

        ApplyRestraints(system, coords, restraintIdxs)
        simulation_initial = Simulation(top.topology, system, m_integrator, platform, properties)
        simulation_initial.context.setPositions(coords.positions)
        return m_integrator, simulation_initial

    def Minimize():
        state = simulation_.context.getState(getEnergy=True)
        print("INITIAL POTENTIAL ENERGY: ", state.getPotentialEnergy())
        print('\nMinimizing local energy...')
        print('Minimizing to 10 kilojoule / (nanometer * mole)')
        simulation_.minimizeEnergy(tolerance=10 * kilojoule / (nanometer * mole), maxIterations=1000)
        print('Minimizing to 1 kilojoule / (nanometer * mole)')
        simulation_.minimizeEnergy(tolerance=1 * kilojoule / (nanometer * mole), maxIterations=1000)
        print('Minimizing to 0.1 kilojoule / (nanometer * mole)')
        simulation_.minimizeEnergy(tolerance=0.1 * kilojoule / (nanometer * mole), maxIterations=1000)
        final_energy = simulation_.context.getState(getEnergy=True)

        print("FINAL POTENTIAL ENERGY: ", final_energy.getPotentialEnergy())

        state = simulation_.context.getState(getPositions=True, getEnergy=True)
        with open('minimized.pdb_', 'w') as output_:
            PDBFile.writeFile(simulation_.topology, state.getPositions(), output_)

    def Thermalize(m_integrator):
        print("\nWarming up the complex with the restraints for 1 ns -> 5e5 steps.")
        steps = 8050
        for t_i in range(5, 311, 5):
            temperature = t_i * kelvin
            m_integrator.setTemperature(temperature)
            simulation_.context.setVelocitiesToTemperature(temperature)
            simulation_.step(steps)
        finalSteps: int = 500000 - 490000
        simulation_.context.setVelocitiesToTemperature(310 * kelvin)
        simulation_.step(finalSteps)

    def EaseRestraints(protocolSteps_):
        easing = int(protocolSteps_ * 0.8)
        remaining_prot = int(protocolSteps_ - easing)
        print(f'\n{easing} steps easing "k" restraint...:')
        easingSteps: int = int(easing / 25)
        for i_T in range(50, 0, -2):
            r_i = (round(i_T * 0.1, 3))
            r_i = max(r_i, 0)
            print('restraint coefficient = ', r_i, "for ", easingSteps, "steps.")
            simulation_.context.setParameter('k', r_i * kilocalories_per_mole / nanometer ** 2)
            simulation_.step(easingSteps)
        print("Completing the remaining protocol steps ", remaining_prot, " with k=0.1")
        simulation_.context.setParameter('k', 0.1)  # safety measure
        simulation_.context.reinitialize(True)
        simulation_.step(remaining_prot)

    # if we're restarting we load the chk and compute the remaining steps
    if parameters['recoverEq'] or parameters['extendEq']:
        integrator, simulation_ = InstantiateSystem()
        if os.path.exists('equilibration_checkpnt.xml'):
            simulation_.loadState('equilibration_checkpnt.xml')

        if os.path.exists('equilibration_checkpnt.chk'):
            simulation_.loadCheckpoint(f'equilibration_checkpnt.chk')
        else:
            raise FileNotFoundError("No checkpoint found to extend the equilibration.")
        step = simulation_.currentStep
        equilibrationSteps: int = int(equilibrationSteps - step)

        if parameters['extendEq']:
            equilibrationSteps = int(((parameters['extendEq'] * 10 ** 6) / parameters['eqIntegrator']) * 0.2)
            protocolSteps = int(((parameters['extendEq'] * 10 ** 6) / parameters['eqIntegrator']) * 0.8)
            e_totalSteps = protocolSteps + equilibrationSteps
            print("Extending equilibration from the last checkpoint. Steps to completion: ", e_totalSteps,
                  f" split in {protocolSteps} protocol steps and {equilibrationSteps} equilibration steps")

    else:  # add easing snippet if it's first time equilibration
        integrator, simulation_ = InstantiateSystem()
        Minimize()
        Thermalize(m_integrator=integrator)
        simulation_.reporters.append(CheckpointReporter('equilibration_checkpnt.chk', parameters['saveFreq']))
        EaseRestraints(protocolSteps)

    print("Setting k=0 and equilibrating for", equilibrationSteps, " steps.")
    simulation_.context.setParameter('k', 0)  # safety measure
    simulation_.context.reinitialize(True)
    simulation_.context.setStepCount(0)
    simulation_.reporters.append(CheckpointReporter('equilibration_checkpnt.chk', parameters['saveFreq']))
    simulation_.reporters.append(
        StateDataReporter('equilibration.std', 1000, step=True, totalSteps=e_totalSteps, speed=True, remainingTime=True, potentialEnergy=True, kineticEnergy=True, temperature=True))
    if os.path.exists('equilibration.xtc') or os.path.exists('equilibration.dcd'):
        append = True
    else:
        append = False
    try:
        simulation_.reporters.append(
            XTCReporter('equilibration.xtc', parameters['saveFreq'], append=append, enforcePeriodicBox=True))
    except:
        simulation_.reporters.append(
            DCDReporter('equilibration.dcd', parameters['saveFreq'], append=append, enforcePeriodicBox=True))
    print("\nEquilibration started.")
    simulation_.step(equilibrationSteps)
    final_state = simulation_.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True,
                                               getParameters=True)
    with open('equilibration_checkpnt.xml', 'w') as output:
        output.write(XmlSerializer.serialize(final_state))
    simulation_.reporters.clear()
    print("\nEquilibration finished successfully.")
