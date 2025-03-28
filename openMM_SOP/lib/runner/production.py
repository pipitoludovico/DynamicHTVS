import os.path
from os import makedirs
from openmm.app import *

from ..builder.mdSettings import GetMDsettings, SetBarostat, ApplyRestraints, RemoveRestraints
from ..FEP.Alchemical import *


def RunProduction(GPU, replica, system, top, parameters, FF, files, coords):
    append, step = False, 0
    productionSteps = int((parameters['productionTime'] * 10 ** 6) / parameters['prIntegrator'])
    friction = parameters['prIntegrator'] / 1000
    integrator, platform, properties = GetMDsettings(GPU, friction)
    makedirs(f"run_{replica}", exist_ok=True)
    simulation_production = Simulation(top.topology, system, integrator, platform, properties)
    # we now check if it's a recover, extension or post-equilibration
    if parameters['recoverPr'] or parameters['extendPr']:  # we load the state if present
        append = True
        if os.path.exists(f'run_{replica}/{replica}_checkpnt.chk'):
            simulation_production.loadCheckpoint(f'run_{replica}/{replica}_checkpnt.chk')
            step = simulation_production.context.getStepCount()
        else:
            raise FileNotFoundError(f"No run_{replica}/{replica}_checkpnt.chk file found from which to restart.")
        if parameters['recoverPr']:  # we count how many steps to completion
            productionSteps: int = int(productionSteps - step)
            print("\nRestarting the production from the last production's checkpoint to completion. Steps to completion: ", step, "/", productionSteps)
        if parameters['extendPr']:
            productionSteps = int((parameters['extendPr'] * 10 ** 6) / parameters['prIntegrator'])
            print("Extending the production by: ", productionSteps, " steps.")
    else:
        if os.path.exists('equilibration_checkpnt.chk'):
            simulation_production.loadCheckpoint('equilibration_checkpnt.chk')

    if parameters['prBarostat']:
        print("Adding barostat to production.")
        SetBarostat(system)
    print("Running for ", productionSteps, "steps")
    simulation_production.context.reinitialize(True)
    simulation_production.context.setStepCount(step)
    simulation_production.reporters.append(CheckpointReporter(f'run_{replica}/{replica}_checkpnt.chk', parameters['saveFreq']))
    try:
        simulation_production.reporters.append(XTCReporter(f"run_{replica}/Traj_{replica}.xtc", parameters['saveFreq'], append=append, enforcePeriodicBox=True))
    except Exception as e:
        print("No XTC writer found. Using a DCD writer\n.", e)
        simulation_production.reporters.append(DCDReporter(f"run_{replica}/Traj_{replica}.dcd", parameters['saveFreq'], append=append, enforcePeriodicBox=True))
    simulation_production.reporters.append(StateDataReporter(f'run_{replica}/Traj_{replica}.std', 1000, step=True, totalSteps=productionSteps, remainingTime=True, potentialEnergy=True, temperature=True, speed=True))

    simulation_production.step(productionSteps)
    final_state = simulation_production.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True, getParameters=True)
    with open(f'run_{replica}/Traj_{replica}.xml', 'w') as output:
        output.write(XmlSerializer.serialize(final_state))

    if parameters['FEP']:
        #
        # prod_sys = simulation_production.context.getSystem()
        # unrestrained = RemoveRestraints(prod_sys)
        # alchemical_system = SetAlchemicalForces(unrestrained, parameters['ligandResname'], top=top, FF=FF, files=files, coords=coords)
        #
        # friction = parameters['prIntegrator'] / 1000
        # integrator, platform, properties = GetMDsettings(GPU, friction)
        #
        # simulation_alchemical = Simulation(top.topology, alchemical_system, integrator, platform, properties)
        #
        # RunFEP(simulation_alchemical, integrator=integrator)
        # CalculateEnergy()
        # In your production.py, modify this section:
        prod_sys = simulation_production.context.getSystem()
        state = simulation_production.context.getState(getPositions=True)
        coords = state.getPositions()
        unrestrained = RemoveRestraints(prod_sys)
        alchemical_system = SetAlchemicalForces(unrestrained, parameters['ligandResname'], top=top, FF=FF, files=files,
                                                coords=coords)

        friction = parameters['prIntegrator'] / 1000
        integrator, platform, properties = GetMDsettings(GPU, friction)

        # Create simulation AND set positions
        simulation_alchemical = Simulation(top.topology, alchemical_system, integrator, platform, properties)
        simulation_alchemical.context.setPositions(coords)  # <-- This was missing

        RunFEP(simulation_alchemical, integrator=integrator)

    simulation_production.reporters.clear()
