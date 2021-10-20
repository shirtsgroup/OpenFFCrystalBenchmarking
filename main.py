# Import Block
import numpy as np
from rdkit import Chem
from openbabel import pybel
from openff.toolkit.topology import Molecule,Topology
from openff.toolkit.utils import RDKitToolkitWrapper
import openmm
from simtk.openmm import app
from openmm import unit
from openmm.app import PDBFile
from openff.toolkit.typing.engines.smirnoff import ForceField
from openmm.vec3 import Vec3
import warnings
warnings.simplefilter("ignore")


# Wrapper and FF setup
# Use RDKit wrapper
rdktkw = RDKitToolkitWrapper()

# Loading setup parameters
forcefield = ForceField('openff-2.0.0.offxml')

#%%

# load pdb with one copy of pdb file
off_mol = Molecule.from_pdb_and_smiles('data/7101899.pdb', "CC1=CN=C(C(=C1OC)C)C[S@@](=O)C2=NC3=C(N2)C=C(C=C3)OC")
# load supercell pdb file (2x2x2) into topology
pdb_file = PDBFile('data/7101899_supercell.pdb')
off_top = Topology.from_openmm(pdb_file.topology, [off_mol])
# Create MD simulation inputs
system = forcefield.create_openmm_system(off_top)
integrator = openmm.VerletIntegrator(1*unit.femtoseconds)
platform = openmm.Platform.getPlatformByName('Reference')

#%%

# create simulation
simulation = openmm.app.Simulation(pdb_file.topology, system, integrator, platform)
# set initial positions from pdbfile
positions = pdb_file.getPositions()
simulation.context.setPositions(positions)

#%%

# set reporters
pdb_reporter = openmm.app.PDBReporter('trajectory.pdb', 1)
state_data_reporter = openmm.app.StateDataReporter(
    "data.csv",
    1,
    step=True,
    potentialEnergy=True,
    temperature=True,
    density=True,
)
simulation.reporters.append(pdb_reporter)
simulation.reporters.append(state_data_reporter)

#%%

simulation.context.setPositions(positions)
simulation.saveState('initial')
orig_potential = simulation.context.getState(getEnergy=True).getPotentialEnergy()
#state_data_reporter.report(simulation, simulation.context.getState())
print('Initial Energy ' + str(orig_potential))
print('Minimizing Energy!')
simulation.minimizeEnergy()
min_state = simulation.context.getState(getEnergy=True, getPositions=True, getForces=True)
min_potential = min_state.getPotentialEnergy()
print('Final Energy = ' + str(min_potential))

#%%

# Report if system uses periodic boundary conditions and forces used
print(simulation.system.usesPeriodicBoundaryConditions())
print(simulation.system.getForces())

#%%

simulation.step(1)


#%%

import pymol
pymol.cmd.load('7101899_supercell.pdb','initial')
pymol.cmd.load('trajectory.pdb','final')
rmsd = pymol.cmd.align('final','initial',cycles=0)
print(rmsd[0])

#%%

# Need block to feed parameters to minimizer, where forces provide 3*n derivatives of energy
# wrt position, +3 more derivatives of energy wrt box vectors
# This block is a work in progress
forces = simulation.context.getState(getForces=True).getForces()
p_box_pos = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()

x0 = np.array([2,0,2,0,0,2]) # a_x, b_x, b_y, c_x, c_y, c_z
def box_energy(sim,x,positions):
    # x will be np array of all inputs
    sim.system.setDefaultPeriodicBoxVectors(Vec3(x[0],0,0),Vec3(x[1],x[2],0),Vec3(x[3],x[4],x[5]))
    print(sim.system.getForces())
    new_integrator = openmm.VerletIntegrator(1*unit.femtoseconds)
    new_sim = openmm.app.Simulation(pdb_file.topology, system, new_integrator, platform)
    new_sim.context.setPositions(positions)
    return new_sim.context.getState(getEnergy=True).getPotentialEnergy()
energy = box_energy(simulation,x0,positions)
print(energy)
#print(forces)
