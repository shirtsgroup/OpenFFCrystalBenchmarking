from rdkit import Chem
from openbabel import pybel
from openff.toolkit.topology import Molecule,Topology
from openff.toolkit.utils import RDKitToolkitWrapper
import openmm
from simtk.openmm import app
from openmm import unit
from openmm.app import PDBFile
from openff.toolkit.typing.engines.smirnoff import ForceField

# Use RDKit wrapper
rdktkw = RDKitToolkitWrapper()

# Loading setup parameters
forcefield = ForceField('openff_unconstrained-1.1.0.offxml')
# load pdb with one copy of pdb file
off_mol = Molecule.from_pdb_and_smiles("data/7101899.pdb", "CC1=CN=C(C(=C1OC)C)C[S@@](=O)C2=NC3=C(N2)C=C(C=C3)OC")
pdb_file = PDBFile('data/7101899_supercell.pdb')
off_top = Topology.from_openmm(pdb_file.topology, [off_mol])
system = forcefield.create_openmm_system(off_top)
integrator = openmm.VerletIntegrator(1*unit.femtoseconds)
platform = openmm.Platform.getPlatformByName('Reference')

# create simulation and minimize energy
simulation = openmm.app.Simulation(pdb_file.topology, system, integrator, platform)
# set initial positions from pdbfile
positions = pdb_file.getPositions()

simulation.context.setPositions(positions)
pdb_reporter = openmm.app.PDBReporter('data/trajectory.pdb', 1)
state_data_reporter = openmm.app.StateDataReporter(
    "data/data.csv",
    1,
    step=True,
    potentialEnergy=True,
    temperature=True,
    density=True,
)

simulation.saveState('initial')

#state_data_reporter.report(simulation, simulation.context.getState())
simulation.reporters.append(pdb_reporter)
simulation.reporters.append(state_data_reporter)
print('Minimizing Energy!')
simulation.minimizeEnergy()
simulation.saveState('after')


print("simulating??")
simulation.step(2)
