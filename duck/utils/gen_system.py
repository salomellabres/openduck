import shutil
from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.topology import Molecule, Topology
from rdkit import Chem
import parmed
from simtk.openmm.app import PDBFile

def generateSMIRNOFFStructureRDK(ligand_file, ff='test_forcefields/smirnoff99Frosst.offxml'):
	"""
	Given an RDKit molecule, create an OpenMM System and use to
	generate a ParmEd structure using the SMIRNOFF forcefield parameters.
	"""
	if ligand_file.endswith('.mol'):
		new_file = ligand_file.replace('.mol', '.sdf')
		shutil.copyfile(ligand_file, new_file)
		ligand_file = new_file
	ligand_off_molecule = Molecule.from_file(ligand_file, allow_undefined_stereo=True)
	force_field = ForceField(ff)
	
	ligand_topology = Topology.from_molecules(ligand_off_molecule)
	# ligand_system = force_field.create_openmm_system(ligand_topology, charge_from_molecules=[ligand_off_molecule])
	ligand_system = force_field.create_openmm_system(ligand_topology)
	ligand_topology = ligand_topology.to_openmm()  # needed for call to parmed
	
	# Read in the coordinates of the ligand from the PDB file
	Chem.MolToPDBFile(Chem.MolFromMolFile(ligand_file, removeHs=False), "ligand.pdb")
	ligand_pdbfile = PDBFile("ligand.pdb")
	# Convert OpenMM System object containing ligand parameters into a ParmEd Structure.
	ligand_structure = parmed.openmm.load_topology(ligand_pdbfile.topology, ligand_system, xyz=ligand_pdbfile.positions)
	return ligand_structure
