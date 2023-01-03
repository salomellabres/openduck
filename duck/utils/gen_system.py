import shutil

from rdkit import Chem
import parmed
from simtk.openmm.app import PDBFile
from openff.toolkit.topology import Molecule, Topology


def generateSMIRNOFFStructureRDK_old(ligand_file, ff='test_forcefields/smirnoff99Frosst.offxml'):
	from openff.typing.engines.smirnoff import ForceField
	"""
	Given an RDKit molecule, create an OpenMM System and use to
	generate a ParmEd structure using the SMIRNOFF forcefield parameters.
	"""
	print('Parametrizing with SMIRNOFF')
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

def generateSMIRNOFFStructureRDK(ligand_file, ff='openff-1.0.0'):
	from openmmforcefields.generators import SMIRNOFFTemplateGenerator
	from simtk.openmm.app import ForceField
	"""
	Given an mol file, create an openMM system and use it to generate a ParmEd structure using GAFF2
	"""
	print('Parametrizing with SMIRNOFF')
	if ligand_file.endswith('.mol'):
		new_file = ligand_file.replace('.mol', '.sdf')
		shutil.copyfile(ligand_file, new_file)
		ligand_file = new_file
	ligand_off_molecule = Molecule.from_file(ligand_file, allow_undefined_stereo=True)
	ligand_topology = Topology.from_molecules(ligand_off_molecule).to_openmm()
	smirnoff = SMIRNOFFTemplateGenerator(molecules=ligand_off_molecule, forcefield=ff)
	force_field = ForceField()
	force_field.registerTemplateGenerator(smirnoff.generator)

	# ligand_system = force_field.create_openmm_system(ligand_topology, charge_from_molecules=[ligand_off_molecule])
	ligand_system = force_field.createSystem(ligand_topology)
	
	# Read in the coordinates of the ligand from the PDB file
	Chem.MolToPDBFile(Chem.MolFromMolFile(ligand_file, removeHs=False), "ligand.pdb")
	ligand_pdbfile = PDBFile("ligand.pdb")
	# Convert OpenMM System object containing ligand parameters into a ParmEd Structure.
	ligand_structure = parmed.openmm.load_topology(ligand_pdbfile.topology, ligand_system, xyz=ligand_pdbfile.positions)
	return ligand_structure

def generateGAFFStructureRDK(ligand_file, ff='gaff-2.11'):
	from openmmforcefields.generators import GAFFTemplateGenerator
	from simtk.openmm.app import ForceField
	"""
	Given an mol file, create an openMM system and use it to generate a ParmEd structure using GAFF2
	"""
	print('Parametrizing with GAFF2')
	if ligand_file.endswith('.mol'):
		new_file = ligand_file.replace('.mol', '.sdf')
		shutil.copyfile(ligand_file, new_file)
		ligand_file = new_file
	ligand_off_molecule = Molecule.from_file(ligand_file, allow_undefined_stereo=True)
	ligand_topology = Topology.from_molecules(ligand_off_molecule).to_openmm()
	gaff = GAFFTemplateGenerator(molecules=ligand_off_molecule, forcefield=ff)
	force_field = ForceField()
	force_field.registerTemplateGenerator(gaff.generator)

	# ligand_system = force_field.create_openmm_system(ligand_topology, charge_from_molecules=[ligand_off_molecule])
	ligand_system = force_field.createSystem(ligand_topology)
	
	# Read in the coordinates of the ligand from the PDB file
	Chem.MolToPDBFile(Chem.MolFromMolFile(ligand_file, removeHs=False), "ligand.pdb")
	ligand_pdbfile = PDBFile("ligand.pdb")
	# Convert OpenMM System object containing ligand parameters into a ParmEd Structure.
	ligand_structure = parmed.openmm.load_topology(ligand_pdbfile.topology, ligand_system, xyz=ligand_pdbfile.positions)
	#gaff.add_moleculeS(ligand_structure)
	return ligand_structure

def generateEspalomaFFStructureRDK(ligand_file, ff='espaloma-0.2.2'):
	from openmmforcefields.generators import EspalomaTemplateGenerator
	from simtk.openmm.app import ForceField
	"""
	Given an mol file, create an openMM system and use it to generate a ParmEd structure using GAFF2
	"""
	print('Parametrizing with Espaloma')
	if ligand_file.endswith('.mol'):
		new_file = ligand_file.replace('.mol', '.sdf')
		shutil.copyfile(ligand_file, new_file)
		ligand_file = new_file
	ligand_off_molecule = Molecule.from_file(ligand_file, allow_undefined_stereo=True)
	ligand_topology = Topology.from_molecules(ligand_off_molecule).to_openmm()
	espaloma = EspalomaTemplateGenerator(molecules=ligand_off_molecule, forcefield=ff)
	force_field = ForceField()
	force_field.registerTemplateGenerator(espaloma.generator)

	# ligand_system = force_field.create_openmm_system(ligand_topology, charge_from_molecules=[ligand_off_molecule])
	ligand_system = force_field.createSystem(ligand_topology)
	
	# Read in the coordinates of the ligand from the PDB file
	Chem.MolToPDBFile(Chem.MolFromMolFile(ligand_file, removeHs=False), "ligand.pdb")
	ligand_pdbfile = PDBFile("ligand.pdb")
	# Convert OpenMM System object containing ligand parameters into a ParmEd Structure.
	ligand_structure = parmed.openmm.load_topology(ligand_pdbfile.topology, ligand_system, xyz=ligand_pdbfile.positions)
	return ligand_structure

def generateFFStructureRDK_from_template(ligand_file, FFtemplate):
	from simtk.openmm.app import ForceField
	"""
	Given an mol file, create an openMM system and use it to generate a ParmEd structure using GAFF2
	"""
	print('Parametrizing with GAFF2')
	if ligand_file.endswith('.mol'):
		new_file = ligand_file.replace('.mol', '.sdf')
		shutil.copyfile(ligand_file, new_file)
		ligand_file = new_file
	ligand_off_molecule = Molecule.from_file(ligand_file, allow_undefined_stereo=True)
	ligand_topology = Topology.from_molecules(ligand_off_molecule).to_openmm()
	gaff = FFtemplate(molecules=ligand_off_molecule)
	force_field = ForceField()
	force_field.registerTemplateGenerator(gaff.generator)

	# ligand_system = force_field.create_openmm_system(ligand_topology, charge_from_molecules=[ligand_off_molecule])
	ligand_system = force_field.createSystem(ligand_topology)
	
	# Read in the coordinates of the ligand from the PDB file
	Chem.MolToPDBFile(Chem.MolFromMolFile(ligand_file, removeHs=False), "ligand.pdb")
	ligand_pdbfile = PDBFile("ligand.pdb")
	# Convert OpenMM System object containing ligand parameters into a ParmEd Structure.
	ligand_structure = parmed.openmm.load_topology(ligand_pdbfile.topology, ligand_system, xyz=ligand_pdbfile.positions)
	return ligand_structure
	
if __name__=='__main__':
	import sys
	r = generateGAFFStructureRDK(sys.argv[1])
