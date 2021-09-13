from rdkit import Chem
from openbabel import pybel
from openff.toolkit.topology import Molecule,Topology
from openff.toolkit.utils import RDKitToolkitWrapper
from simtk import unit
import openmm
from openmm.app import PDBFile
from openff.toolkit.typing.engines.smirnoff import ForceField

# Use RDKit wrapper
rdktkw = RDKitToolkitWrapper()

# Loading setup parameters
forcefield = ForceField('openff_unconstrained-1.1.0.offxml')
off_mol = Molecule.from_pdb_and_smiles("notebooks/2107460.pdb", "O=C(N)N1c2ccccc2C=Cc2ccccc12")
pdb_file = PDBFile('notebooks/2107460.pdb')
off_top = Topology.from_openmm(pdb_file.topology, [off_mol])
system = forcefield.create_openmm_system(off_top)
integrator = openmm.VerletIntegrator(1*unit.femtoseconds)
platform = openmm.Platform.getPlatformByName('Reference')

# create simulation and minimize energy
#simulation = openmm.app.Simulation(pdb_file.topology, system, integrator, platform)
# set initial positions from pdbfile
positions = pdb_file.getPositions()

#simulation.context.setPositions(positions)

#simulation.saveState('before')
#simulation.minimizeEnergy()
#simulation.saveState('after')
