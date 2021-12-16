# Import Block
import os
import numpy as np
from rdkit import Chem
from openff.toolkit.topology import Molecule,Topology
from openff.toolkit.utils import RDKitToolkitWrapper
import openmm
from simtk.openmm import app
from openmm import unit
from openmm.app import PDBFile
from openff.toolkit.typing.engines.smirnoff import ForceField
import mdtraj
import warnings
import logging
import pandas as pd
from scipy.optimize import minimize
import csv

# Functions to minimize energy with respect to box vectors
# forces provide 3*n derivatives of energy
# wrt position, +6 more derivatives of energy wrt box vectors
def box_energy(x,*args):
    # x is np array of 3*n(atoms) positional coordinates and 6 coordinates of 3 triclinical box vectors
    # Return the energy of the system (in kJ/mol)
    # *args will have the simulation context and n (number of particles)
    context = args[0]
    n = args[1]
    # Build position array
    positions_arr = np.empty([n,3])
    for i in range(int(n)):
        positions_arr[i][0] = x[i*3]
        positions_arr[i][1] = x[i*3+1]
        positions_arr[i][2] = x[i*3+2]
    # Build periodic box vectors
    a = np.array([x[n*3],0,0])
    b = np.array([x[n*3+1],x[n*3+2],0])
    c = np.array([x[n*3+3],x[n*3+4],x[n*3+5]])
    # Set Context with positions and periodic boundary conditions
    context.setPositions(positions_arr)
    context.setPeriodicBoxVectors(a,b,c)
    # Return Energy
    energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return energy

def jacobian(x,*args):
    #must return a n*3 + 6 size vector with derivative of energy with respect to each input parameter
    #for positions, return given forces
    context = args[0]
    n = args[1]

    energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    forces = -1*context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/(unit.nano*unit.meter))
    jac = forces.flatten()
    # Positions
    positions_arr = np.empty([n,3])
    for i in range(int(n)):
        positions_arr[i][:] = x[i*3:i*3+3]
    context.setPositions(positions_arr)
    # box vectors
    a = np.array([x[n*3],0,0])
    b = np.array([x[n*3+1],x[n*3+2],0])
    c = np.array([x[n*3+3],x[n*3+4],x[n*3+5]])
    epsilon = 1e-5
    box_vectors = np.array([a,b,c])
    frac_positions = np.matmul(np.linalg.inv(box_vectors),positions_arr.T).T

    # a_x
    temp_a_x = x[n*3] + epsilon
    temp_box_vectors = np.array([np.array([temp_a_x,0,0]),b,c]) # Replace appropriate entry in box vectors
    new_positions = np.matmul(temp_box_vectors,frac_positions.T).T # Transform frac coordiantes with new box vectors
    temp_x = np.concatenate((new_positions.flatten(),[temp_a_x],x[n*3+1:n*3+6])) # build x array to feed to box_energy function with new box vector
    new_energy = box_energy(temp_x,context,n) #compute new energy
    grad_temp = (new_energy-energy)/epsilon # compute gradient
    jac = np.append(jac,grad_temp) # append to jac

    # b_x
    temp_b_x = x[n*3+1] + epsilon
    temp_box_vectors = np.array([a,np.array([temp_b_x,x[n*3+2],0]),c])
    new_positions = np.matmul(temp_box_vectors,frac_positions.T).T
    temp_x = np.concatenate((new_positions.flatten(),[x[n*3]],[temp_b_x],x[n*3+2:n*3+6]))
    new_energy = box_energy(temp_x,context,n)
    grad_temp = (new_energy-energy)/epsilon
    jac = np.append(jac,grad_temp)

    # b_y
    temp_b_y = x[n*3+2] + epsilon
    temp_box_vectors = np.array([a,np.array([x[n*3+1],temp_b_y,0]),c])
    new_positions = np.matmul(temp_box_vectors,frac_positions.T).T
    temp_x = np.concatenate((new_positions.flatten(),x[n*3:n*3+2],[temp_b_y],x[n*3+3:n*3+6]))
    new_energy = box_energy(temp_x,context,n)
    grad_temp = (new_energy-energy)/epsilon
    jac = np.append(jac,grad_temp)

    #c_x
    temp_c_x = x[n*3+3] + epsilon
    temp_box_vectors = np.array([a,b,np.array([temp_c_x,x[n*3+4],x[n*3+5]])])
    new_positions = np.matmul(temp_box_vectors,frac_positions.T).T
    temp_x = np.concatenate((new_positions.flatten(),x[n*3:n*3+3],[temp_c_x],x[n*3+4:n*3+6]))
    new_energy = box_energy(temp_x,context,n)
    grad_temp = (new_energy-energy)/epsilon
    jac = np.append(jac,grad_temp)

    #c_y
    temp_c_y = x[n*3+4] + epsilon
    temp_box_vectors = np.array([a,b,np.array([x[n*3+3],temp_c_y,x[n*3+5]])])
    new_positions = np.matmul(temp_box_vectors,frac_positions.T).T
    temp_x = np.concatenate((new_positions.flatten(),x[n*3:n*3+4],[temp_c_y],[x[n*3+5]]))
    new_energy = box_energy(temp_x,context,n)
    grad_temp = (new_energy-energy)/epsilon
    jac = np.append(jac,grad_temp)

    #c_z
    temp_c_z = x[n*3+5] + epsilon
    temp_box_vectors = np.array([a,b,np.array([x[n*3+3],x[n*3+4],temp_c_z])])
    new_positions = np.matmul(temp_box_vectors,frac_positions.T).T
    temp_x = np.concatenate((new_positions.flatten(),x[n*3:n*3+5],[temp_c_z]))
    new_energy = box_energy(temp_x,context,n)
    grad_temp = (new_energy-energy)/epsilon
    jac = np.append(jac,grad_temp)

    return jac

# Warning and Logger Setup
warnings.simplefilter("ignore")
logging.basicConfig(filename='errors.log', filemode='w')

# Wrapper and FF setup
# Use RDKit wrapper
rdktkw = RDKitToolkitWrapper()

# Loading setup parameters
forcefield = ForceField('openff-1.0.0.offxml')
# Load smiles csv file as pandas
smiles = pd.read_csv('allcod.smi', names=['SMILES', 'COD ID'], sep='\t')
# Initialize text file for rmsd values to be recorded
with open('example_data/rmsd_values.txt','w') as f:
    f.write('COD ID\tRMSD\n')
# Initialize empty array to store data from sucessful minimizations
data = []
for pdb in os.listdir('example_data/PDB'):
    # load pdb with one copy of pdb file
    cod_id = pdb.split('.')[0]
    print(cod_id)
    try:
        try:
            # get smiles from all_smilies
            smiles_string = smiles.loc[smiles['COD ID'] == int(cod_id)]['SMILES'].values[0]
            off_mol = Molecule.from_pdb_and_smiles('example_data/PDB/' + pdb, smiles_string)
        except Exception as e:
            logging.error('PDB/SMILES error with ID %s' % cod_id)
            logging.error(e)
            continue
        try:
            # load supercell pdb file into topology
            pdb_file = PDBFile('example_data/PDB_supercell/' + cod_id + '_supercell.pdb')
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
        pdb_reporter = openmm.app.PDBReporter('example_data/minimized_PDB_supercell/' + cod_id + '.pdb', 1)
        simulation.reporters.append(pdb_reporter)
        # set positions
        simulation.context.setPositions(positions)
        # save state and print initial PE
        simulation.saveState('example_data/initial_states/' + cod_id + '_initial.xml')
        orig_potential = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        if orig_potential.value_in_unit(unit.kilojoule_per_mole) > 1e24:
            #Skip if the initial energy evaluation is infeasibly high
            continue
        print('Initial Energy ' + str(orig_potential))
        # Minimize Energy and save final state
        print('Minimizing Energy!')
        simulation.minimizeEnergy(maxIterations=100000) # Perform initial minimization with openMM minimizer
        mid_potential = simulation.context.getState(getEnergy=True).getPotentialEnergy()
        # # build x array to feed to minimizer
        box_vectors = pdb_file.topology.getPeriodicBoxVectors()
        A = np.array(box_vectors.value_in_unit(unit.nano * unit.meter))
        numpy_positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
        x = np.append(numpy_positions.flatten(), [A[0][0], A[1][0], A[1][1], A[2][0], A[2][1], A[2][2]])
        n = len(numpy_positions)
        # run minimizer
        result = minimize(box_energy, x, (simulation.context, n), method='L-BFGS-B', jac=jacobian, options={'maxiter': 100})
        print(result)
        x_new = result.x
        # # update simulation with minimized positions and box vectors
        positions_arr = np.empty([n, 3])
        for i in range(int(n)):
           positions_arr[i][:] = x_new[i * 3:i * 3 + 3]
        simulation.context.setPositions(positions_arr)
        a = np.array([x_new[n * 3], 0, 0])
        b = np.array([x_new[n * 3 + 1], x_new[n * 3 + 2], 0])
        c = np.array([x_new[n * 3 + 3], x_new[n * 3 + 4], x_new[n * 3 + 5]])
        simulation.context.setPeriodicBoxVectors(a, b, c)

        min_state = simulation.context.getState(getEnergy=True, getPositions=True, getForces=True)
        min_potential = min_state.getPotentialEnergy()
        simulation.saveState('example_data/final_states/' + cod_id + '_final.xml')
        print('Final Energy = ' + str(min_potential))

        simulation.step(1)

        initial = mdtraj.load_pdb('example_data/PDB_supercell/' + cod_id + '_supercell.pdb')
        final = mdtraj.load_pdb('example_data/minimized_PDB_supercell/' + cod_id + '.pdb')
        rmsd = mdtraj.rmsd(initial, final)
        data.append(
            {
                'COD ID': cod_id,
                'RMSD': rmsd,
                'Original Energy': orig_potential,
                'Minimized Energy (OpenMM)': mid_potential,
                'Minimized Energy (Cell Minimization)': min_potential,
                'Original Box Vectors': A,
                'Minimized Box Vectors': np.array([a, b, c])
            })
        with open('example_data/rmsd_values.txt','a') as f:
            f.write('%s\t%s\n' % (cod_id, rmsd[0]))
    except Exception as e:
        logging.error('Generic error with ID %s' % cod_id)
        logging.error(e)
        continue

try:
    d = pd.DataFrame(data)
    d.to_csv('example_data/minimization_results.csv')
    d.to_pickle('example_data/minimization_results.pkl')
except Exception as e:
    logging.error('Error with save of results data')



