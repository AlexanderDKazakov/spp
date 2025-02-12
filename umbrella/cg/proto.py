__author__ = "Alexander D. Kazakov (aleksandr[dot]kazakov_at_uv[dot]es)"
__version__ = "0.1.0"

from config import CONFIG, OUTPUT_PATH, GROMACS_TOP, MOTIONS
from config import ENERGY_DECOMPOSITION_FNAME

import openmm as mm
import openmm.app as app
import martini_openmm as martini
import numpy as np
from tqdm import tqdm

if __name__ == "__main__":

    gro = app.GromacsGroFile('input.gro')
    top = martini.MartiniTopFile(
        "input.top",
        periodicBoxVectors=gro.getPeriodicBoxVectors(),
        defines={},
        epsilon_r=15.0,
        includeDir=GROMACS_TOP,
    )

    # Creating system
    system = top.create_system(nonbonded_cutoff = CONFIG['nonBondedCutoff'])
    integrator = mm.LangevinMiddleIntegrator(CONFIG['temperature'], CONFIG['frictionCoeff'], CONFIG['step_size'])

    # Energy contribution: this must be done before you create a Simulation.
    for i, f in enumerate(system.getForces()): f.setForceGroup(i)

    simulation = app.Simulation(top.topology, system, integrator)
    simulation.context.setPositions(gro.positions)

    # minimization
    simulation.minimizeEnergy()
    state = simulation.context.getState(getPositions=True, getEnergy=True,)
    with open(OUTPUT_PATH/"minimized.pdb", "w") as output: app.PDBFile.writeFile(simulation.topology, state.getPositions(), output)

    # starting value
    desired_value = CONFIG['motion_type_params'][CONFIG['motion_type']]['initial_variable_value']
    for idxs in CONFIG['motion_type_params'][CONFIG['motion_type']]['atom_indices']:

        if CONFIG['motion_type'] == MOTIONS.TRANSLATION:
            for idx1, idx2 in [idxs]:
                cv = mm.CustomBondForce('r')
                cv.addBond(idx1, idx2)

                pullingForce = mm.CustomCVForce(CONFIG['motion_type_params'][CONFIG['motion_type']]['cv_formula'])
                pullingForce.addGlobalParameter("fc_pull", CONFIG['motion_type_params'][CONFIG['motion_type']]['fc_pull'])
                pullingForce.addGlobalParameter(CONFIG['motion_type_params'][CONFIG['motion_type']]['variable'], desired_value)
                pullingForce.addCollectiveVariable("cv", cv)
                system.addForce(pullingForce)

        elif CONFIG['motion_type'] == MOTIONS.ROTATION:
            for idx1, idx2, idx3, idx4 in [idxs]:
                cv = mm.CustomTorsionForce("theta")
                cv.addTorsion(idx1, idx2, idx3, idx4)
                pullingForce = mm.CustomCVForce(CONFIG['motion_type_params'][CONFIG['motion_type']]['cv_formula'])
                pullingForce.addGlobalParameter("fc_pull", CONFIG['motion_type_params'][CONFIG['motion_type']]['fc_pull'])
                pullingForce.addGlobalParameter(CONFIG['motion_type_params'][CONFIG['motion_type']]['variable'], desired_value)
                pullingForce.addCollectiveVariable("cv", cv)
                system.addForce(pullingForce)

    simulation.context.reinitialize(preserveState=True)

    window_coords = []
    window_index = 0
    windows = CONFIG['motion_type_params'][CONFIG['motion_type']]['windows']


    with open(OUTPUT_PATH/ENERGY_DECOMPOSITION_FNAME, "w") as decompose:
        decompose.write("# IDX WIN ") # header
        for idx in range(2):
            for i, f in enumerate(system.getForces()):
                state = simulation.context.getState(getEnergy=True, groups={i})
                # print(f.getName(), state.getPotentialEnergy())
                if idx == 0: decompose.write(f.getName() + f"[{str(state.getPotentialEnergy()).split()[1]}]" + " ") # name with units
                else:
                    if i == 0: decompose.write(f"\n{str(simulation.currentStep)} -1 ")
                    decompose.write(str(state.getPotentialEnergy()._value) + " ")

    # SMD pulling loop
    for i in range(CONFIG['total_steps']//CONFIG['increment_steps']):
        simulation.step(CONFIG['increment_steps'])
        current_cv_value = pullingForce.getCollectiveVariableValues(simulation.context)

        current_desired_ratio = (1 - abs(desired_value/current_cv_value[0])._value) * 100 # %
        if (i*CONFIG['increment_steps']) % 5000 == 0:
            print(f"desired = {desired_value}| current value = {current_cv_value} || {current_desired_ratio:.2f}%")

        # increment the location of the CV based on the pulling velocity
        desired_value += CONFIG['motion_type_params'][CONFIG['motion_type']]['velocity_pulling'] * CONFIG['step_size'] * CONFIG['increment_steps']
        simulation.context.setParameter(CONFIG['motion_type_params'][CONFIG['motion_type']]['variable'], desired_value)

        # check if we should save this config as a window starting structure
        if (window_index < len(windows) and current_desired_ratio < 5 and current_cv_value >= windows[window_index]):
            print(f"[Window:{window_index}] desired value = {desired_value}| current value = {current_cv_value[0]}")
            window_coords.append(simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions())
            window_index += 1

        if window_index == len(windows): break

    # save the window structures
    for i, coords in enumerate(window_coords):
        with open(OUTPUT_PATH / f'window_{i}.pdb', 'w') as outfile:
            app.PDBFile.writeFile(simulation.topology, coords, outfile)

    # run windows
    for idx, desired_value in enumerate(windows):
        # load starting PDB file
        with open(OUTPUT_PATH / f'window_{idx}.pdb') as infile: pdb = app.PDBFile(infile)

        # reusing the existing simulation
        simulation.context.setPositions(pdb.positions)
        simulation.context.setParameter(CONFIG['motion_type_params'][CONFIG['motion_type']]['variable'], desired_value)
        simulation.context.setVelocitiesToTemperature(CONFIG['temperature'])

        simulation.reporters.append(app.PDBReporter(OUTPUT_PATH/f"output_{idx}.pdb", CONFIG['report_steps'], enforcePeriodicBox=False))
        simulation.reporters.append(app.StateDataReporter(str(OUTPUT_PATH/f"output_{idx}.dat"), CONFIG['report_steps'], step=True, potentialEnergy=True, temperature=True, totalEnergy=True))

        # run short equilibration with new positions and desired_value
        print(f"Window [{idx}] | Short equilibration:{CONFIG['windows_short_equilib_steps']*CONFIG['step_size']}")
        simulation.step(CONFIG['windows_short_equilib_steps'])

        # run the data collection
        cv_values = []
        print(f"Window [{idx}] | Total run:{CONFIG['total_steps']*CONFIG['step_size']}")
        for i in tqdm(range(CONFIG['total_steps']//CONFIG['cv_record_steps'])):
            try:
                simulation.step(CONFIG['cv_record_steps'])
            except ValueError as e:
                print("Hmmm strange step! --> ", e)
                continue

            # getting the current cv
            current_cv_value = pullingForce.getCollectiveVariableValues(simulation.context)
            cv_values.append([i, current_cv_value[0]])
            # energy decomposition
            with open(OUTPUT_PATH/ENERGY_DECOMPOSITION_FNAME, "a") as decompose:
                for i, f in enumerate(system.getForces()):
                    state = simulation.context.getState(getEnergy=True, groups={i})
                    # print(f.getName(), state.getPotentialEnergy())
                    if i == 0: decompose.write(f"\n{str(simulation.currentStep)} {str(idx)} ")
                    decompose.write(str(state.getPotentialEnergy()._value) + " ")

        np.savetxt(OUTPUT_PATH/f"cv_values_window_{idx}.txt", np.array(cv_values))
        print(f"Completed window [{idx}]")

    print("Done")

