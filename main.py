# Import Block
import os
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
import mdtraj
import warnings
import logging
import pandas as pd

# Warning and Logger Setup
warnings.simplefilter("ignore")
logging.basicConfig(filename='errors.log', filemode='w')

# Wrapper and FF setup
# Use RDKit wrapper
rdktkw = RDKitToolkitWrapper()

# Loading setup parameters
forcefield = ForceField('openff-2.0.0.offxml')
# Load smiles csv file as pandas
smiles = pd.read_csv('allcod.smi', names=['SMILES', 'COD ID'], sep='\t')
with open('data/rmsd_values.txt','w') as f:
    f.write('COD ID\tRMSD\n')

for pdb in os.listdir('data/PDB'):
    try:
        # load pdb with one copy of pdb file
        cod_id = pdb.split('.')[0]
        print(cod_id)
        try:
            # get smiles from all_smilies
            smiles_string = smiles.loc[smiles['COD ID'] == int(cod_id)]['SMILES'].values[0]
            off_mol = Molecule.from_pdb_and_smiles('data/PDB/' + pdb, smiles_string)
        except Exception as e:
            logging.error('PDB/SMILES error with ID %s' % cod_id)
            logging.error(e)
            continue
        try:
            # load supercell pdb file (2x2x2) into topology
            pdb_file = PDBFile('data/PDB_supercell/' + cod_id + '_supercell.pdb')
            off_top = Topology.from_openmm(pdb_file.topology, [off_mol])
        except Exception as e:
            logging.error('Topology error with ID %s' % cod_id)
            logging.error(e)
            continue
        # Create MD simulation inputs
        system = forcefield.create_openmm_system(off_top)
        integrator = openmm.VerletIntegrator(1*unit.femtoseconds)
        platform = openmm.Platform.getPlatformByName('Reference')

        try:
            # create simulation, catch errors
            simulation = openmm.app.Simulation(pdb_file.topology, system, integrator, platform)
        except Exception as e:
            logging.error('Simulation Build error with ID %s' % cod_id)
            logging.error(e)
            continue
        # set initial positions from pdbfile
        positions = pdb_file.getPositions()
        simulation.context.setPositions(positions)

        # set reporters
        pdb_reporter = openmm.app.PDBReporter('data/trajectories/' + cod_id + '.pdb', 1)
        simulation.reporters.append(pdb_reporter)
        # set positions
        simulation.context.setPositions(positions)
        # save state and print initial PE
        simulation.saveState('data/initial_states/' + cod_id + '_initial.xml')
        orig_potential = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        if orig_potential.value_in_unit(unit.kilojoule_per_mole) > 1e24:
            continue
        print('Initial Energy ' + str(orig_potential))
        # Minimize Energy and save final state
        print('Minimizing Energy!')
        simulation.minimizeEnergy(maxIterations=100000)
        min_state = simulation.context.getState(getEnergy=True, getPositions=True, getForces=True)
        min_potential = min_state.getPotentialEnergy()
        simulation.saveState('data/final_states/' + cod_id + '_final.xml')
        print('Final Energy = ' + str(min_potential))

        simulation.step(1)

        initial = mdtraj.load_pdb('data/PDB_supercell/' + cod_id + '_supercell.pdb')
        final = mdtraj.load_pdb('data/trajectories/' + cod_id + '.pdb')
        rmsd = mdtraj.rmsd(initial, final)
        with open('data/rmsd_values.txt','a') as f:
            f.write('%s\t%s\n' % (cod_id, rmsd[0]))
    except Exception as e:
        logging.error('Generic Error')
        logging.error(e)
        continue
