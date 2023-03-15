import simtk.openmm as mm
import simtk.openmm.app as app
import simtk.unit as u
import sys, os
import pickle
from duck.utils import duck_stuff
from duck.utils import cal_ints


def perform_md(
    checkpoint_in_file,
    checkpoint_out_file,
    csv_out_file,
    pdb_out_file,
    dcd_out_file,
    force_constant_ligand=1.0,
    md_len=1.0,
    force_constant_chunk=0.1,
    gpu_id=0,
):
    '''Performs molecular dynamics (MD) simulation of a ligand-protein complex using OpenMM.

Arguments:
- checkpoint_in_file (str): Path to the input checkpoint file.
- checkpoint_out_file (str): Path to the output checkpoint file.
- csv_out_file (str): Path to the CSV output file.
- pdb_out_file (str): Path to the PDB output file.
- dcd_out_file (str): Path to the DCD output file.
- force_constant_ligand (float, optional): Force constant of the harmonic restraint on the ligand position (in kJ/mol/nm^2). Default is 1.0.
- md_len (float, optional): Length of the MD simulation (in nanoseconds). Default is 1.0.
- force_constant_chunk (float, optional): Force constant of the harmonic restraint on the chunk position (in kJ/mol/nm^2). Default is 0.1.
- gpu_id (int, optional): ID of the GPU device to use, or None for CPU. Default is 0.

If the output checkpoint file already exists, the function returns without performing the simulation.

The input checkpoint file should be generated by a previous call to this function, or by some other means of initializing the simulation state.

The function expects a file called 'complex_system.pickle' in the current working directory, which should contain a serialized OpenMM System object representing the ligand-protein complex.

The function calculates the heavy atoms in the complex, applies a harmonic restraint to their positions, and applies a restraint to the distance between the ligand and a selected interaction site.

The MD simulation is performed using a Langevin integrator with a timestep of 2 femtoseconds. The simulation state is saved to the output checkpoint file, and the final coordinates are saved to the output PDB file.

The function also generates a CSV file containing various state data of the simulation, and a DCD file containing the trajectory of the simulation.

Raises:
- FileNotFoundError: If the 'complex_system.pickle' file is not found.
- ImportError: If OpenMM or one of its dependencies cannot be imported.
- RuntimeError: If an OpenMM-related error occurs during the simulation.
'''

    if os.path.isfile(checkpoint_out_file):
        print(f'{checkpoint_out_file} is already calculated, skipping' )
        return
    #print("loading pickle")
    pickle_in = open("complex_system.pickle", "rb")
    pkl = pickle.load(pickle_in)
    combined_pmd = pkl[0]
    # print(dir(combined_pmd)) # why print the dir of the combined system from the pickle?
    key_interaction = pkl[1:]
    pickle_in.close()
    MD_len = md_len * u.nanosecond
    sim_steps = round(MD_len / (0.002 * u.picosecond))
    # Platform definition
    
    platformProperties = {}
    if gpu_id != None:
        platform = mm.Platform_getPlatformByName("CUDA")
        platformProperties["CudaPrecision"] = "double"
        platformProperties["DeviceIndex"] = str(gpu_id)
    else:
        platform = mm.Platform_getPlatformByName("CPU")
    platformProperties["DeterministicForces"] = 'true'
    # Get indexes of heavy atoms in chunk
    Chunk_Heavy_Atoms = duck_stuff.getHeavyAtomsInSystem(combined_pmd)
    # Setting System
    system = combined_pmd.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=9 * u.angstrom,
        constraints=app.HBonds,
        hydrogenMass=None,
    )
    # Apply force on all havy atoms of chunk and apply restraint for the ligand-chunk distance
    duck_stuff.applyHarmonicPositionalRestraints(
        system, force_constant_chunk, combined_pmd.positions, Chunk_Heavy_Atoms
    )
    duck_stuff.applyLigandChunkRestraint(
        system,
        force_constant_ligand,
        10.0,
        2 * u.angstrom,
        3 * u.angstrom,
        4 * u.angstrom,
        key_interaction,
    )
    # Integrator
    integrator = mm.LangevinIntegrator(
        300 * u.kelvin, 4 / u.picosecond, 0.002 * u.picosecond
    )
    # Setting Simulation object and loading the checkpoint
    simulation = app.Simulation(
        combined_pmd.topology, system, integrator, platform, platformProperties
    )
    simulation.loadCheckpoint(checkpoint_in_file)
    # Simulation reporters
    simulation.reporters.append(
        app.StateDataReporter(
            csv_out_file,
            2000,
            step=True,
            time=True,
            totalEnergy=True,
            kineticEnergy=True,
            potentialEnergy=True,
            temperature=True,
            density=True,
            progress=True,
            totalSteps=sim_steps,
            speed=True,
        )
    )
    simulation.reporters.append(app.DCDReporter(dcd_out_file, 2000))

    # Production

    simulation.step(sim_steps)
    # Save state in checkpoint file and save coordinates in PDB file
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    positions = state.getPositions()
    app.PDBFile.writeFile(simulation.topology, positions, open(pdb_out_file, "w"))
    simulation.saveCheckpoint(checkpoint_out_file)


if __name__ == "__main__":
    if len(sys.argv) != 7:
        sys.exit("Usage 03_md.py in.chk out.chk out.csv out.pdb")
    checkpoint_in_file = sys.argv[1]
    checkpoint_out_file = sys.argv[2]
    csv_out_file = sys.argv[3]
    pdb_out_file = sys.argv[4]
    perform_md(checkpoint_in_file, checkpoint_out_file, csv_out_file, pdb_out_file)
