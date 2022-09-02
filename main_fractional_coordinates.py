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
import pdb
from openff.toolkit.typing.engines.smirnoff import ForceField
import mdtraj
import warnings
import logging
import pandas as pd
from scipy.optimize import minimize, basinhopping
import csv
from pyxtal.operations import angle, create_matrix
from pyxtal.constants import deg, rad, ltype_keywords
import matplotlib.pyplot as plt


# Functions to minimize energy with respect to box parameter
# forces provide 3*n derivatives of energy
# wrt frac position, +6 more derivatives of energy wrt box parameter
# def convert_fractional_coordinate():


def matrix2para(matrix, radians=False):
    """
    Given a 3x3 matrix representing a unit cell, outputs a list of lattice
    parameters.
    Args:
        matrix: a 3x3 array or list, where the first, second, and third rows
            represent the a, b, and c vectors respectively
        radians: if True, outputs angles in radians. If False, outputs in
            degrees
    Returns:
        a 1x6 list of lattice parameters [a, b, c, alpha, beta, gamma]. a, b,
        and c are the length of the lattice vectos, and alpha, beta, and gamma
        are the angles between these vectors (in radians by default)
    """
    cell_para = np.zeros(6)
    # a
    cell_para[0] = np.linalg.norm(matrix[0])
    # b
    cell_para[1] = np.linalg.norm(matrix[1])
    # c
    cell_para[2] = np.linalg.norm(matrix[2])
    # alpha
    cell_para[3] = angle(matrix[1], matrix[2])
    # beta
    cell_para[4] = angle(matrix[0], matrix[2])
    # gamma
    cell_para[5] = angle(matrix[0], matrix[1])

    if not radians:
        # convert radians to degrees
        deg = 180.0 / np.pi
        cell_para[3] *= deg
        cell_para[4] *= deg
        cell_para[5] *= deg
    return cell_para

def para2matrix(cell_para, radians=False, format="lower"):
    """
    Given a set of lattic parameters, generates a matrix representing the
    lattice vectors
    Args:
        cell_para: a 1x6 list of lattice parameters [a, b, c, alpha, beta,
            gamma]. a, b, and c are the length of the lattice vectos, and
            alpha, beta, and gamma are the angles between these vectors. Can
            be generated by matrix2para
        radians: if True, lattice parameters should be in radians. If False,
            lattice angles should be in degrees
        format: a string ('lower', 'symmetric', or 'upper') for the type of
            matrix to be output
    Returns:
        a 3x3 matrix representing the unit cell. By default (format='lower'),
        the a vector is aligined along the x-axis, and the b vector is in the
        y-z plane
    """
    a = cell_para[0]
    b = cell_para[1]
    c = cell_para[2]
    alpha = cell_para[3]
    beta = cell_para[4]
    gamma = cell_para[5]
    if radians is not True:
        alpha *= rad
        beta *= rad
        gamma *= rad
    cos_alpha = np.cos(alpha)
    cos_beta = np.cos(beta)
    cos_gamma = np.cos(gamma)
    sin_gamma = np.sin(gamma)
    sin_alpha = np.sin(alpha)
    matrix = np.zeros([3, 3])
    if format == "lower":
        # Generate a lower-diagonal matrix
        c1 = c * cos_beta
        c2 = (c * (cos_alpha - (cos_beta * cos_gamma))) / sin_gamma
        matrix[0][0] = a
        matrix[1][0] = b * cos_gamma
        matrix[1][1] = b * sin_gamma
        matrix[2][0] = c1
        matrix[2][1] = c2
        matrix[2][2] = np.sqrt(c ** 2 - c1 ** 2 - c2 ** 2)
    elif format == "symmetric":
        # TODO: allow generation of symmetric matrices
        pass
    elif format == "upper":
        # Generate an upper-diagonal matrix
        a3 = a * cos_beta
        a2 = (a * (cos_gamma - (cos_beta * cos_alpha))) / sin_alpha
        matrix[2][2] = c
        matrix[1][2] = b * cos_alpha
        matrix[1][1] = b * sin_alpha
        matrix[0][2] = a3
        matrix[0][1] = a2
        tmp = a ** 2 - a3 ** 2 - a2 ** 2
        if tmp > 0:
            matrix[0][0] = np.sqrt(a ** 2 - a3 ** 2 - a2 ** 2)
        else:
            return None
        #pass
    return matrix




def box_energy(x,*args):
    # x is np array of 3*n(atoms) fractional position and 6 triclinical box parameters
    # Return the energy of the system (in kJ/mol)
    # *args will have the simulation context and n (number of particles)
    context = args[0]
    n = args[1]
    # Build frac position array
    frac_positions_arr = np.empty([n,3])
    for i in range(int(n)):
        frac_positions_arr[i][0] = x[i*3]
        frac_positions_arr[i][1] = x[i*3+1]
        frac_positions_arr[i][2] = x[i*3+2]
    # Build 6 triclinical box parameters
    box_parameter = [x[n*3], x[n*3+1], x[n*3+2], x[n*3+3], x[n*3+4], x[n*3+5]]
    # convert box parameter to box vector
    box_vectors = para2matrix(box_parameter, radians=False, format="lower")
    # Build periodic box vectors
    a = np.array([box_vectors[0][0],0,0])
    b = np.array([box_vectors[1][0],box_vectors[1][1],0])
    c = np.array([box_vectors[2][0],box_vectors[2][1],box_vectors[2][2]])
    # (TEST) print to monitor the result
    #print("-------------------")
    #print("periodic box vectors")
    #print(a)
    #print(b)
    #print(c)
    #(TEST) convert fractional position to atom position
    positions_arr = np.matmul(box_vectors, frac_positions_arr.T).T
    # Set Context with positions and periodic boundary conditions
    context.setPositions(positions_arr)
    context.setPeriodicBoxVectors(a,b,c)
    # Return Energy
    energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return energy


# (plot) list for storing data to plot
grad_temp_func = []
grad_temp_func0 = []
grad_temp_func1 = []
grad_temp_func2 = []
grad_temp_func3 = []
grad_temp_func4 = []
grad_temp_func5 = []
value_0 = []
value_1 = []
value_2 = []
value_3 = []
value_4 = []
value_5 = []
reduced_positions_gradient_x0 = []
reduced_positions_gradient_y0 = []
reduced_positions_gradient_z0 = []
reduced_positions_gradient_x100 = []
reduced_positions_gradient_y100 = []
reduced_positions_gradient_z100 = []
reduced_positions_gradient_x200 = []
reduced_positions_gradient_y200 = []
reduced_positions_gradient_z200 = []


def jacobian(x,*args):
    # x is np array of 3*n(atoms) fractional position and 6 triclinical box parameters
    # must return a n*3 + 6 box parameters with derivative of energy with respect to each input parameter
    # for positions, return given forces
    context = args[0]
    n = args[1]
    # energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    forces = -1*context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole/(unit.nano*unit.meter))
    jac = forces.flatten()

    # Positions
    frac_positions_arr = np.empty([n,3])
    for i in range(int(n)):
        frac_positions_arr[i][:] = x[i*3:i*3+3]

    # convert box parameter to box vector for atom fractional position
    box_parameter = [x[n * 3], x[n * 3 + 1], x[n * 3 + 2], x[n * 3 + 3], x[n * 3 + 4], x[n * 3 + 5]]
    box_vectors = para2matrix(box_parameter, radians=False, format="lower")

    #(TEST) convert fractional position to atom position
    positions_arr = np.matmul(box_vectors, frac_positions_arr.T).T
    context.setPositions(positions_arr)


    # (plot) Store reduced positions gradient vs iteration
    reduced_positions_gradient = np.matmul(box_vectors, forces.T).T
    reduced_positions_gradient_x0.append(reduced_positions_gradient[0][0])
    reduced_positions_gradient_y0.append(reduced_positions_gradient[0][1])
    reduced_positions_gradient_z0.append(reduced_positions_gradient[0][2])
    reduced_positions_gradient_x100.append(reduced_positions_gradient[100][0])
    reduced_positions_gradient_y100.append(reduced_positions_gradient[100][1])
    reduced_positions_gradient_z100.append(reduced_positions_gradient[100][2])
    reduced_positions_gradient_x200.append(reduced_positions_gradient[200][0])
    reduced_positions_gradient_y200.append(reduced_positions_gradient[200][1])
    reduced_positions_gradient_z200.append(reduced_positions_gradient[200][2])

    # A, B, C, alpha, beta, gamma
    for i in range(len(box_parameter)):
        # use different epsilon for A, B, C and alpha, beta, gamma
        epsilon = 1e-6
        temp_U = x[n * 3 + i] + epsilon * x[n * 3 + i]
        temp_L = x[n * 3 + i] - epsilon * x[n * 3 + i]
        dstep = temp_U - temp_L

        temp_box_parameters = np.array([x[n * 3], x[n * 3 + 1], x[n * 3 + 2], x[n * 3 + 3], x[n * 3 + 4], x[n * 3 + 5]])
        temp_box_parameters[i] = temp_U # Replace appropriate entry in box parameter
        temp_box_vectors = para2matrix(temp_box_parameters, radians=False, format="lower") # convert box parameter to box vector for atom fractional position
        # (TEST) Transform frac coordiantes with new box vectors
        new_positions = np.matmul(temp_box_vectors, frac_positions_arr.T).T
        # (TEST) convert atom position to fractional position
        new_frac_positions_arr = np.matmul(np.linalg.inv(temp_box_vectors), new_positions.T).T
        temp_x = np.concatenate((new_frac_positions_arr.flatten(), temp_box_parameters.flatten()))
        energy_U = box_energy(temp_x, context, n)  # compute new energy

        temp_box_parameters = np.array([x[n * 3], x[n * 3 + 1], x[n * 3 + 2], x[n * 3 + 3], x[n * 3 + 4], x[n * 3 + 5]])
        temp_box_parameters[i] = temp_L # Replace appropriate entry in box parameter
        temp_box_vectors = para2matrix(temp_box_parameters, radians=False, format="lower") # convert box parameter to box vector for atom fractional position
        # (TEST) Transform frac coordiantes with new box vectors
        new_positions = np.matmul(temp_box_vectors, frac_positions_arr.T).T
        # (TEST) convert atom position to fractional position
        new_frac_positions_arr = np.matmul(np.linalg.inv(temp_box_vectors), new_positions.T).T
        temp_x = np.concatenate((new_frac_positions_arr.flatten(), temp_box_parameters.flatten()))
        energy_L = box_energy(temp_x, context, n)  # compute new energy
        grad_temp = (energy_U - energy_L) / dstep # compute gradient
        jac = np.append(jac, grad_temp)  # append to jac

        # (plot) store gradient/value_temp vs iteration
        if i == 0:
            grad_temp_func0.append(grad_temp)
            value_0.append(x[n * 3 + i])
        elif i == 1:
            grad_temp_func1.append(grad_temp)
            value_1.append(x[n * 3 + i])
        elif i == 2:
            grad_temp_func2.append(grad_temp)
            value_2.append(x[n * 3 + i])
        elif i == 3:
            grad_temp_func3.append(grad_temp)
            value_3.append(x[n * 3 + i])
        elif i == 4:
            grad_temp_func4.append(grad_temp)
            value_4.append(x[n * 3 + i])
        elif i == 5:
            grad_temp_func5.append(grad_temp)
            value_5.append(x[n * 3 + i])

    return jac




# Warning and Logger Setup
warnings.simplefilter("ignore")
logging.basicConfig(filename='errors.log', filemode='w')

# Wrapper and FF setup
# Use RDKit wrapper
rdktkw = RDKitToolkitWrapper()

# Loading setup parameters

forcefield=ForceField("openff_unconstrained-2.0.0.offxml")

# Load smiles csv file as pandas
smiles = pd.read_csv('allcod.smi', names=['SMILES', 'COD ID'], sep='\t')
# Initialize text file for rmsd values to be recorded
with open('data/rmsd_values.txt','w') as f:
    f.write('COD ID\tRMSD\n')
# Initialize empty array to store data from sucessful minimizations
data = []
for pdb in os.listdir('data/PDB'):
    # load pdb with one copy of pdb file
    cod_id = pdb.split('.')[0]
    print(cod_id)
#    try:
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
        # replace cod_id to
        pdb_file = PDBFile('data/PDB_supercell/' + cod_id + '_supercell.pdb')
        off_top = Topology.from_openmm(pdb_file.topology, [off_mol])
    except Exception as e:
        logging.error('Topology error with ID %s' % cod_id)
        logging.error(e)
        continue
    # Create MD simulation inputs
    system = forcefield.create_openmm_system(off_top)
    integrator = openmm.VerletIntegrator(1 * unit.femtoseconds)
    platform = openmm.Platform.getPlatformByName('CUDA')

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
    pdb_reporter = openmm.app.PDBReporter('data/minimized_PDB_supercell/' + cod_id + '.pdb', 1)
    simulation.reporters.append(pdb_reporter)
    # set positions
    simulation.context.setPositions(positions)
    # save state and print initial PE
    simulation.saveState('data/initial_states/' + cod_id + '_initial.xml')
    orig_potential = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    if orig_potential.value_in_unit(unit.kilojoule_per_mole) > 1e24:
        # Skip if the initial energy evaluation is infeasibly high
        continue
    print('Initial Energy ' + str(orig_potential))
    # Minimize Energy and save final state
    print('Minimizing Energy!')
    simulation.minimizeEnergy(maxIterations=100000)  # Perform initial minimization with openMM minimizer
    mid_potential = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    # build x array to feed to minimizer
    box_vectors = pdb_file.topology.getPeriodicBoxVectors()
    A = np.array(box_vectors.value_in_unit(unit.nano * unit.meter))
    numpy_positions = simulation.context.getState(getPositions=True).getPositions(asNumpy=True)
    print("Initial box vector\n", A)
    # convert box vector to box parameter
    box_parameter = matrix2para(A, radians=False)
    print("Initial box parameter\n", box_parameter)
    box_vectors = para2matrix(box_parameter, radians=False, format="lower")
    # convert atom position to fractional position
    frac_positions_arr = np.matmul(np.linalg.inv(box_vectors), numpy_positions.T).T
    x = np.append(frac_positions_arr.flatten(), [box_parameter[0], box_parameter[1], box_parameter[2], box_parameter[3], box_parameter[4], box_parameter[5]])
    n = len(frac_positions_arr)


    # constraint for minimizer
    def constraint1(x):
        box_parameter = [x[n * 3], x[n * 3 + 1], x[n * 3 + 2], x[n * 3 + 3], x[n * 3 + 4], x[n * 3 + 5]]
        box_vector = para2matrix(box_parameter, radians=False, format="lower")
        box_vector.flatten()
        value = box_vector[0][0] - 2 * abs(box_vector[1][0])
        return value
    def constraint2(x):
        box_parameter = [x[n * 3], x[n * 3 + 1], x[n * 3 + 2], x[n * 3 + 3], x[n * 3 + 4], x[n * 3 + 5]]
        box_vector = para2matrix(box_parameter, radians=False, format="lower")
        value = box_vector[0][0] - 2 * abs(box_vector[2][0])
        return value
    def constraint3(x):
        box_parameter = [x[n * 3], x[n * 3 + 1], x[n * 3 + 2], x[n * 3 + 3], x[n * 3 + 4], x[n * 3 + 5]]
        box_vector = para2matrix(box_parameter, radians=False, format="lower")
        value = box_vector[1][1] - 2 * abs(box_vector[2][1])
        return value
    con1 = {'type': 'ineq', 'fun': constraint1}
    con2 = {'type': 'ineq', 'fun': constraint2}
    con3 = {'type': 'ineq', 'fun': constraint3}
    cons = ([con1, con2, con3])

    # run minimizer
    try:
        # try L-BFGS-B minimization first
        method = 'L-BFGS-B'
        result = minimize(box_energy, x, (simulation.context, n), method='L-BFGS-B', jac=jacobian, options={'disp': 1, 'maxiter': 1000})
    except:
        # try trust-constr minimization first
        method = 'trust-constr'
        result = minimize(box_energy, x, (simulation.context, n), method='trust-constr', constraints=cons, jac=jacobian, options={'disp': 1, 'xtol': 1e-08, 'gtol': 1e-08, 'maxiter': 1000})




    x_new = result.x
    print("result", x_new)
    # update simulation with minimized positions and box vectors
    frac_positions_arr = np.empty([n, 3])
    for i in range(int(n)):
        frac_positions_arr[i][:] = x_new[i * 3:i * 3 + 3]


    # convert result to box vector
    box_parameter_new = [x_new[n * 3], x_new[n * 3 + 1], x_new[n * 3 + 2], x_new[n * 3 + 3], x_new[n * 3 + 4], x_new[n * 3 + 5]]
    box_vectors_new = para2matrix(box_parameter_new, radians=False, format="lower")
    a = np.array([box_vectors_new[0][0], 0, 0])
    b = np.array([box_vectors_new[1][0], box_vectors_new[1][1], 0])
    c = np.array([box_vectors_new[2][0], box_vectors_new[2][1], box_vectors_new[2][2]])
    simulation.context.setPeriodicBoxVectors(a, b, c)

    # convert fractional position to atom position
    positions_arr = np.matmul(box_vectors_new, frac_positions_arr.T).T
    simulation.context.setPositions(positions_arr)

    min_state = simulation.context.getState(getEnergy=True, getPositions=True, getForces=True)
    min_potential = min_state.getPotentialEnergy()
    simulation.saveState('data/final_states/' + cod_id + '_final.xml')
    print('Final Energy = ' + str(min_potential))

    simulation.step(1)

    initial = mdtraj.load_pdb('data/PDB_supercell/' + cod_id + '_supercell.pdb')
    final = mdtraj.load_pdb('data/minimized_PDB_supercell/' + cod_id + '.pdb')
    rmsd = mdtraj.rmsd(initial, final)
    data.append(
        {
            'COD ID': cod_id,
            'RMSD': rmsd,
            'method': method,
            'Original Energy': orig_potential,
            'Minimized Energy (OpenMM)': mid_potential,
            'Minimized Energy (Cell Minimization)': min_potential,
            'Original Box Vectors': A,
            'Minimized Box Vectors': np.array([a, b, c]),
            'Original Box parameters': box_parameter,
            'Minimized Box parameter': box_parameter_new,
            'Box parameter - Last five gradient': " ",
            'Gradient A': grad_temp_func0[-5:-1],
            'Gradient B': grad_temp_func1[-5:-1],
            'Gradient C': grad_temp_func2[-5:-1],
            'Gradient alpha': grad_temp_func3[-5:-1],
            'Gradient beta': grad_temp_func4[-5:-1],
            'Gradient gamma': grad_temp_func5[-5:-1],
            'Atom position - Last five gradient': " ",
            'Gradient reduced positions (x0)': reduced_positions_gradient_x0[-5:-1],
            'Gradient reduced positions (y0)': reduced_positions_gradient_y0[-5:-1],
            'Gradient reduced positions (z0)': reduced_positions_gradient_z0[-5:-1],
            'Gradient reduced positions (x100)': reduced_positions_gradient_x100[-5:-1],
            'Gradient reduced positions (y100)': reduced_positions_gradient_y100[-5:-1],
            'Gradient reduced positions (z100)': reduced_positions_gradient_z100[-5:-1],
            'Gradient reduced positions (x200)': reduced_positions_gradient_x200[-5:-1],
            'Gradient reduced positions (y200)': reduced_positions_gradient_y200[-5:-1],
            'Gradient reduced positions (z200)': reduced_positions_gradient_z200[-5:-1],

        })
    with open('data/rmsd_values.txt', 'a') as f:
        f.write('%s\t%s\n' % (cod_id, rmsd[0]))
#    except Exception as e:
#        logging.error('Generic error with ID %s' % cod_id)
#        logging.error(e)
#        continue

# (plot) gradient + parameter vs iteration
    plt.switch_backend('agg')
    plt.ylim(-10000, 10000)
    plt.ylabel('gradient', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('gradient(A) vs iteration', fontsize=15)
    plt.plot(grad_temp_func0)
    plt.tight_layout()
    plt.savefig('./figure/A_gradient.jpg')
    plt.clf()
    plt.ylabel('parameter', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('parameter(A) vs iteration', fontsize=15)
    plt.plot(value_0)
    plt.tight_layout()
    plt.savefig('./figure/A.jpg')
    plt.clf()


    plt.ylim(-10000, 10000)
    plt.ylabel('gradient', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('gradient(B) vs iteration', fontsize=15)
    plt.plot(grad_temp_func1)
    plt.tight_layout()
    plt.savefig('./figure/B_gradient.jpg')
    plt.clf()
    plt.ylabel('parameter', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('parameter(B) vs iteration', fontsize=15)
    plt.plot(value_1)
    plt.tight_layout()
    plt.savefig('./figure/B.jpg')
    plt.clf()

    plt.ylim(-10000, 10000)
    plt.ylabel('gradient', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('gradient(C) vs iteration', fontsize=15)
    plt.plot(grad_temp_func2)
    plt.tight_layout()
    plt.savefig('./figure/C_gradient.jpg')
    plt.clf()
    plt.ylabel('parameter', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('parameter(C) vs iteration', fontsize=15)
    plt.plot(value_2)
    plt.tight_layout()
    plt.savefig('./figure/C.jpg')
    plt.clf()


    plt.ylim(-10000, 10000)
    plt.ylabel('gradient', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('gradient(alpha) vs iteration', fontsize=15)
    plt.plot(grad_temp_func3)
    plt.tight_layout()
    plt.savefig('./figure/alpha_gradient.jpg')
    plt.clf()
    plt.ylabel('parameter', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('parameter(alpha) vs iteration', fontsize=15)
    plt.plot(value_3)
    plt.tight_layout()
    plt.savefig('./figure/alpha.jpg')
    plt.clf()

    plt.tight_layout()
    plt.ylim(-10000, 10000)
    plt.ylabel('gradient', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('gradient(beta) vs iteration', fontsize=15)
    plt.plot(grad_temp_func4)
    plt.savefig('./figure/beta_gradient.jpg')
    plt.clf()
    plt.ylabel('parameter', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('parameter(beta) vs iteration', fontsize=15)
    plt.plot(value_4)
    plt.savefig('./figure/beta.jpg')
    plt.clf()

    plt.tight_layout()
    plt.ylim(-10000, 10000)
    plt.ylabel('gradient', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('gradient(gamma) vs iteration', fontsize=15)
    plt.plot(grad_temp_func5)
    plt.tight_layout()
    plt.savefig('./figure/gamma_gradient.jpg')
    plt.clf()
    plt.ylabel('parameter', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('parameter(gamma) vs iteration', fontsize=15)
    plt.plot(value_5)
    plt.tight_layout()
    plt.savefig('./figure/gamma.jpg')
    plt.clf()


# (plot)  reduced positions gradient vs iteration

    plt.tight_layout()
    plt.ylim(-10000, 10000)
    plt.ylabel('gradient', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('gradient(gamma) vs iteration', fontsize=15)
    plt.plot(grad_temp_func5)
    plt.tight_layout()
    plt.savefig('./figure/gamma_gradient.jpg')
    plt.clf()
    plt.ylabel('parameter', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('parameter(gamma) vs iteration', fontsize=15)
    plt.plot(value_5)
    plt.tight_layout()
    plt.savefig('./figure/gamma.jpg')
    plt.clf()

    plt.tight_layout()
    plt.ylabel('reduced positions gradient x0', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('reduced positions gradient(x0) vs iteration', fontsize=15)
    plt.plot(reduced_positions_gradient_x0)
    plt.tight_layout()
    plt.savefig('./figure/x0.jpg')
    plt.clf()

    plt.tight_layout()
    plt.ylabel('reduced positions gradient y0', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('reduced positions gradient(y0) vs iteration', fontsize=15)
    plt.plot(reduced_positions_gradient_y0)
    plt.tight_layout()
    plt.savefig('./figure/y0.jpg')
    plt.clf()

    plt.tight_layout()
    plt.ylabel('reduced positions gradient z0', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('reduced positions gradient(z0) vs iteration', fontsize=15)
    plt.plot(reduced_positions_gradient_z0)
    plt.tight_layout()
    plt.savefig('./figure/z0.jpg')
    plt.clf()

    plt.tight_layout()
    plt.ylabel('reduced positions gradient x100', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('reduced positions gradient(x100) vs iteration', fontsize=15)
    plt.plot(reduced_positions_gradient_x100)
    plt.tight_layout()
    plt.savefig('./figure/x100.jpg')
    plt.clf()

    plt.tight_layout()
    plt.ylabel('reduced positions gradient y100', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('reduced positions gradient(y100) vs iteration', fontsize=15)
    plt.plot(reduced_positions_gradient_y100)
    plt.tight_layout()
    plt.savefig('./figure/y100.jpg')
    plt.clf()

    plt.tight_layout()
    plt.ylabel('reduced positions gradient z100', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('reduced positions gradient(z100) vs iteration', fontsize=15)
    plt.plot(reduced_positions_gradient_z100)
    plt.tight_layout()
    plt.savefig('./figure/z100.jpg')
    plt.clf()

    plt.tight_layout()
    plt.ylabel('reduced positions gradient x200', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('reduced positions gradient(x200) vs iteration', fontsize=15)
    plt.plot(reduced_positions_gradient_x200)
    plt.tight_layout()
    plt.savefig('./figure/x200.jpg')
    plt.clf()

    plt.tight_layout()
    plt.ylabel('reduced positions gradient y200', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('reduced positions gradient(y200) vs iteration', fontsize=15)
    plt.plot(reduced_positions_gradient_y200)
    plt.tight_layout()
    plt.savefig('./figure/y200.jpg')
    plt.clf()

    plt.tight_layout()
    plt.ylabel('reduced positions gradient z200', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('reduced positions gradient(z200) vs iteration', fontsize=15)
    plt.plot(reduced_positions_gradient_z200)
    plt.tight_layout()
    plt.savefig('./figure/z200.jpg')
    plt.clf()


try:
    d = pd.DataFrame(data)
    d.to_csv('data/minimization_results.csv')
    d.to_pickle('data/minimization_results.pkl')
except Exception as e:
    logging.error('Error with save of results data')


'''
loss_func = []
def callback(x):
    fobj = box_energy(x,simulation.context, n)
    loss_func.append(fobj)
    iteration = [i for i in range(len(loss_func))]
    plt.plot(iteration, loss_func)
    plt.ylabel('loss_func', fontsize=14)
    plt.xlabel('iteration', fontsize=14)
    plt.title('loss_func vs iteration', fontsize=15)
    plt.show()
    callback=callback,
'''

'''
bounds = []
# atom position no limit
for i in range(n * 3):
    bounds.append((None, None))
# box parameter > 0
for i in range(len(box_parameter)):
    bounds.append((0, None))
'''


'''
cons = [{'type': 'ineq', 'fun': lambda x: x[n * 3] - 2 * x[n * 3 + 1]},
        {'type': 'ineq', 'fun': lambda x: x[n * 3] - 2 * x[n * 3 + 3]},
        {'type': 'ineq', 'fun': lambda x: x[n * 3 + 2] - 2 * x[n * 3 + 4]}]
constraints = cons
'''

'''
# plot grad_temp vs iteration
plt.ylim(-10000, 10000)
plt.plot(grad_temp_func, color='b')
plt.ylabel('grad_temp_func', fontsize=12)
plt.xlabel('iteration', fontsize=12)
plt.title('grad_temp vs iteration', fontsize=15)
plt.tight_layout()
plt.show()
'''