# Import Block
import os
import numpy as np
from rdkit import Chem
from openff.toolkit.topology import Molecule, Topology
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
import copy
import re


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
    return matrix


def box_energy(x, *args):
    # x is np array of 3*n(atoms) fractional position and 6 triclinical box parameters
    # Return the energy of the system (in kJ/mol)
    # *args will have the simulation context and n (number of particles)
    context = args[0]
    n = args[1]
    # Build frac position array
    frac_positions_arr = np.empty([n, 3])
    for i in range(int(n)):
        frac_positions_arr[i][0] = x[i * 3]
        frac_positions_arr[i][1] = x[i * 3 + 1]
        frac_positions_arr[i][2] = x[i * 3 + 2]
    # Build 6 triclinical box parameters
    box_parameter = [x[n * 3], x[n * 3 + 1], x[n * 3 + 2], x[n * 3 + 3], x[n * 3 + 4], x[n * 3 + 5]]
    # convert box parameter to box vector
    box_vectors = para2matrix(box_parameter, radians=False, format="lower")
    # Build periodic box vectors
    a = np.array([box_vectors[0][0], 0, 0])
    b = np.array([box_vectors[1][0], box_vectors[1][1], 0])
    c = np.array([box_vectors[2][0], box_vectors[2][1], box_vectors[2][2]])
    # convert fractional position to atom position
    positions_arr = np.matmul(box_vectors, frac_positions_arr.T).T
    # Set Context with positions and periodic boundary conditions
    context.setPositions(positions_arr)
    context.setPeriodicBoxVectors(a, b, c)
    # Return Energy
    energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return energy


def jacobian(x, *args):
    # x is np array of 3*n(atoms) fractional position and 6 triclinical box parameters
    # must return a n*3 + 6 box parameters with derivative of energy with respect to each input parameter
    # for positions, return given forces
    context = args[0]
    n = args[1]
    # energy = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    forces = -1 * context.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / (unit.nano * unit.meter))
    jac = forces.flatten()

    # Positions
    frac_positions_arr = np.empty([n, 3])
    for i in range(int(n)):
        frac_positions_arr[i][:] = x[i * 3:i * 3 + 3]

    # convert box parameter to box vector for atom fractional position
    box_parameter = [x[n * 3], x[n * 3 + 1], x[n * 3 + 2], x[n * 3 + 3], x[n * 3 + 4], x[n * 3 + 5]]
    box_vectors = para2matrix(box_parameter, radians=False, format="lower")

    # convert fractional position to atom position
    positions_arr = np.matmul(box_vectors, frac_positions_arr.T).T
    context.setPositions(positions_arr)

    # A, B, C, alpha, beta, gamma
    for i in range(len(box_parameter)):
        # use different epsilon for A, B, C and alpha, beta, gamma
        epsilon = 1e-6
        temp_U = x[n * 3 + i] + epsilon * x[n * 3 + i]
        temp_L = x[n * 3 + i] - epsilon * x[n * 3 + i]
        dstep = temp_U - temp_L

        temp_box_parameters = np.array([x[n * 3], x[n * 3 + 1], x[n * 3 + 2], x[n * 3 + 3], x[n * 3 + 4], x[n * 3 + 5]])
        temp_box_parameters[i] = temp_U  # Replace appropriate entry in box parameter
        temp_box_vectors = para2matrix(temp_box_parameters, radians=False,
                                       format="lower")  # convert box parameter to box vector for atom fractional position
        # Transform frac coordiantes with new box vectors
        new_positions = np.matmul(temp_box_vectors, frac_positions_arr.T).T
        # convert atom position to fractional position
        new_frac_positions_arr = np.matmul(np.linalg.inv(temp_box_vectors), new_positions.T).T
        temp_x = np.concatenate((new_frac_positions_arr.flatten(), temp_box_parameters.flatten()))
        energy_U = box_energy(temp_x, context, n)  # compute new energy

        temp_box_parameters = np.array([x[n * 3], x[n * 3 + 1], x[n * 3 + 2], x[n * 3 + 3], x[n * 3 + 4], x[n * 3 + 5]])
        temp_box_parameters[i] = temp_L  # Replace appropriate entry in box parameter
        temp_box_vectors = para2matrix(temp_box_parameters, radians=False,
                                       format="lower")  # convert box parameter to box vector for atom fractional position
        # Transform frac coordiantes with new box vectors
        new_positions = np.matmul(temp_box_vectors, frac_positions_arr.T).T
        # convert atom position to fractional position
        new_frac_positions_arr = np.matmul(np.linalg.inv(temp_box_vectors), new_positions.T).T
        temp_x = np.concatenate((new_frac_positions_arr.flatten(), temp_box_parameters.flatten()))
        energy_L = box_energy(temp_x, context, n)  # compute new energy
        grad_temp = (energy_U - energy_L) / dstep  # compute gradient
        jac = np.append(jac, grad_temp)  # append to jac
    return jac


# Tobias RMSD
def get_rmsd(pdb_path_1, pdb_path_2):
    """Calculate Remove Periodicity RMSD

    keyword arguments:
    pdb_path_1(list)   -- the list of atom position before minimization
    pdb_path_2(list)   -- the list of atom position after minimization

    return:
    rmsd(float)        -- RMSD
    """
    if not os.path.exists(pdb_path_1):
        return None
    if not os.path.exists(pdb_path_2):
        return None
    pos1 = app.PDBFile(pdb_path_1).getPositions(asNumpy=True)
    pos2 = app.PDBFile(pdb_path_2).getPositions(asNumpy=True)
    topology = app.PDBFile(pdb_path_1).getTopology()
    non_H_idxs = list()
    for atom_idx, atom in enumerate(topology.atoms()):
        if atom.element.atomic_number != 1:
            non_H_idxs.append(atom_idx)

    diff = np.linalg.norm(
        pos1[non_H_idxs] - pos2[non_H_idxs],
        axis=1
    )
    rmsd = np.sqrt(
        np.mean(
            diff ** 2
        )
    )
    return rmsd


def rmv_periodicity(pos1Array, pos2Array):
    """Calculate Remove Periodicity RMSD

    keyword arguments:
    pos1Array(Array)   -- the array of atom position before minimization
    pos2Array(Array)   -- the array of atom position after minimization

    return
    pos2Array_r(Array)   -- the array that had been removed Periodicity
    """

    pos2Array_r = pos2Array
    for i in range(len(pos1Array)):
        for j in range(3):
            diff = abs((pos1Array[i][j] - pos2Array_r[i][j]))
            if diff > 0.5 * box_parameter_new[j]:
                if pos2Array_r[i][j] > pos1Array[i][j]:
                    tmp = pos2Array_r[i][j] - box_parameter_new[j]
                elif pos2Array_r[i][j] < pos1Array[i][j]:
                    tmp = pos2Array_r[i][j] + box_parameter_new[j]
                else:
                    tmp = pos2Array_r[i][j]
                pos2Array[i][j] = tmp

    return pos2Array_r


# Remove Periodicity RMSD
def get_rmsd_r(pdb_path_1, pdb_path_2):
    """Calculate Remove Periodicity RMSD

    keyword arguments:
    pdb_path_1(list)   -- the list of atom position before minimization
    pdb_path_2(list)   -- the list of atom position after minimization

    return:
    rmsd(float)        -- RMSD that remove periodicity
    """

    if not os.path.exists(pdb_path_1):
        return None
    if not os.path.exists(pdb_path_2):
        return None

    # pos1Array (exp) 2D Array - 1. atom 2. x y z
    # pos2Array (simulation) 2D Array - 1. atom 2. x y z
    pos1Array = np.array(app.PDBFile(pdb_path_1).getPositions(asNumpy=True))
    pos2Array = np.array(app.PDBFile(pdb_path_2).getPositions(asNumpy=True))

    # Record non hydrogen and remove hydrogen for RMSD calculation
    topology = app.PDBFile(pdb_path_1).getTopology()
    non_H_idxs = list()
    for atom_idx, atom in enumerate(topology.atoms()):
        if atom.element.atomic_number != 1:
            non_H_idxs.append(atom_idx)

    # remove periodicity
    pos2Array_r = rmv_periodicity(pos1Array, pos2Array)

    diff = np.linalg.norm(
        pos1Array[non_H_idxs] - pos2Array_r[non_H_idxs],
        axis=1
    )

    rmsd = np.sqrt(
        np.mean(
            diff ** 2
        )
    )

    return rmsd


# RMSD 20
def get_rmsd_n(pdb_path_1, pdb_path_2, option):
    """Calculate Remove Periodicity RMSD

    keyword arguments:
    pdb_path_1(list)   -- the list of atom position before minimization
    pdb_path_2(list)   -- the list of atom position after minimization

    return:
    option = 0
    rmsd(float)        -- RMSD
    option = r
    rmsd(float)        -- RMSD that remove periodicity
    option = int
    rmsd(list)         -- RMSD for option molecule

    """

    if not os.path.exists(pdb_path_1):
        return None
    if not os.path.exists(pdb_path_2):
        return None

    # pos1Array (exp) 2D Array - 1. atom 2. x y z
    # pos2Array (simulation) 2D Array - 1. atom 2. x y z
    pos1Array = np.array(app.PDBFile(pdb_path_1).getPositions(asNumpy=True))
    pos2Array = np.array(app.PDBFile(pdb_path_2).getPositions(asNumpy=True))

    # Record non hydrogen and remove hydrogen for RMSD calculation
    topology = app.PDBFile(pdb_path_1).getTopology()
    non_H_idxs = list()
    for atom_idx, atom in enumerate(topology.atoms()):
        if atom.element.atomic_number != 1:
            non_H_idxs.append(atom_idx)

    # Tobias RMSD
    if option == "o":
        diff = np.linalg.norm(
            pos1Array[non_H_idxs] - pos2Array[non_H_idxs],
            axis=1
        )
        rmsd = np.sqrt(
            np.mean(
                diff ** 2
            )
        )

    # Remove Periodicity RMSD
    elif option == "r":

        # remove periodicity
        pos2Array_r = rmv_periodicity(pos1Array, pos2Array)

        diff = np.linalg.norm(
            pos1Array[non_H_idxs] - pos2Array_r[non_H_idxs],
            axis=1
        )

        rmsd = np.sqrt(
            np.mean(
                diff ** 2
            )
        )



    # RMSD_20 Calculation
    elif isinstance(option, int):

        n = option
        rmsd = np.zeros(int(number_of_molecule_in_supercell))

        # remove periodicity
        pos2Array_r = rmv_periodicity(pos1Array, pos2Array)

        # Remove periodicity for calculation of distance. From first molecule, second, ...
        for target in range(int(number_of_molecule_in_supercell)):
            pos1Array_mol = pos1Array.copy()
            pos2Array_mol = pos2Array_r.copy()
            pos1Array_mol = pos1Array_mol.reshape(int(number_of_molecule_in_supercell), number_of_atom, 3)
            pos2Array_mol = pos2Array_mol.reshape(int(number_of_molecule_in_supercell), number_of_atom, 3)

            for i in range(int(number_of_molecule_in_supercell)):
                for j in range(number_of_atom):
                    for k in range(3):
                        diff = abs((pos2Array_mol[target][j][k] - pos2Array_mol[i][j][k]))
                        if diff > 0.5 * box_parameter_new[k]:
                            if pos2Array_mol[target][j][k] < pos2Array_mol[i][j][k]:
                                tmp1 = pos1Array_mol[i][j][k] - box_parameter[k]
                                tmp2 = pos2Array_mol[i][j][k] - box_parameter_new[k]
                            elif pos2Array_mol[target][j][k] > pos2Array_mol[i][j][k]:
                                tmp1 = pos1Array_mol[i][j][k] + box_parameter[k]
                            else:
                                tmp1 = pos1Array_mol[i][j][k]
                                tmp2 = pos2Array_mol[i][j][k]
                            pos1Array_mol[i][j][k] = tmp1
                            pos2Array_mol[i][j][k] = tmp2

            # Calculate average molecule position
            mol_pos = np.mean(pos2Array_mol, axis=1)

            # Distance calculation
            mol_dis = np.zeros(int(number_of_molecule_in_supercell))
            for i in range(int(number_of_molecule_in_supercell)):
                mol_dis_diff = np.array(mol_pos[i]) - np.array(mol_pos[target])
                distance = mol_dis_diff[0] ** 2 + mol_dis_diff[1] ** 2 + mol_dis_diff[2] ** 2
                mol_dis[i] = distance

            # sort and find the index of n the closest molecule
            close_idx = np.argsort(mol_dis)[:n].tolist()

            # Calculate the diff for n closest molecule
            n_pos_diff = np.zeros((int(number_of_molecule_in_supercell), int(number_of_atom), 3))
            for i in range(len(pos1Array_mol)):
                if i in close_idx:
                    n_pos_diff[i] = pos1Array_mol[i] - pos2Array_mol[i]

            # reshape format
            n_pos_diff = n_pos_diff.flatten().reshape(int(number_of_molecule_in_supercell) * number_of_atom, 3)

            diff = np.linalg.norm(
                n_pos_diff[non_H_idxs],
                axis=1
            )

            rmsd_target = np.sqrt(
                np.mean(
                    diff ** 2
                )
            )

            rmsd[target] = rmsd_target

    return rmsd


# Warning and Logger Setup
warnings.simplefilter("ignore")
logging.basicConfig(filename='errors.log', filemode='w')

# Wrapper and FF setup
# Use RDKit wrapper
rdktkw = RDKitToolkitWrapper()

# Loading setup parameters
forcefield = ForceField('openff_unconstrained-2.0.0.offxml')
# Load smiles csv file as pandas
smiles = pd.read_csv('allcod.smi', names=['SMILES', 'COD ID'], sep='\t')
# Initialize text file for rmsd values to be recorded
with open('data/rmsd_values.txt', 'w') as f:
    f.write('COD ID\tRMSD\n')
# Initialize empty array to store data from sucessful minimizations
data = []
for pdb in os.listdir('data/PDB/2100147'):
    # load pdb with one copy of pdb file
    cod_id = pdb.split('.')[0]
    print(cod_id)

    # try:
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
        pdb_file = PDBFile('data/PDB_supercell/' + "2100147" + '_supercell.pdb')
        off_top = Topology.from_openmm(pdb_file.topology, [off_mol])
    except Exception as e:
        logging.error('Topology error with ID %s' % cod_id)
        logging.error(e)
        continue

    # Create MD simulation inputs and GPU setting
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
    x = np.append(frac_positions_arr.flatten(),
                  [box_parameter[0], box_parameter[1], box_parameter[2], box_parameter[3], box_parameter[4],
                   box_parameter[5]])
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

    # Extract number of molecules in unitcell
    cif_file = open('data/CIF/' + cod_id + '.cif', 'r')
    for line in cif_file.readlines():
        if line.startswith("_cell_formula_units_Z"):
            number_of_molecule = int(re.findall(r"\d+", line)[0])

    # Extract number of atoms in molecule
    cif_file = open('data/CIF/' + cod_id + '.cif', 'r')
    for line in cif_file.readlines():
        if line.startswith("_chemical_formula_sum"):
            chemical_formula = line.split(' ')[-10:]
            number_of_atom = 0
            for i in range(len(chemical_formula)):
                if len(chemical_formula[i]) == 0:
                    number_of_atom += 0
                elif len(chemical_formula[i]) == 1:
                    number_of_atom += 1
                else:
                    try:
                        number_of_atom += int(re.findall(r"\d+", chemical_formula[i])[0])
                    except:
                        number_of_atom += 1

    # Extract number of atoms in supercell
    f_supercell = open('data/PDB_supercell/' + cod_id + '_supercell.pdb', 'r')
    number_of_atom_in_supercell = 0
    for line in f_supercell.readlines():
        if line.startswith("HETATM"):
            number_of_atom_in_supercell += 1

    # Calculate number of molecules in supercell
    number_of_molecule_in_supercell = number_of_atom_in_supercell / number_of_atom

    # run minimizer
    try:
        # try L-BFGS-B minimization first
        method = 'L-BFGS-B'
        result = minimize(box_energy, x, (simulation.context, n), method='L-BFGS-B', jac=jacobian,
                          options={'disp': 1, 'maxiter': 1})
    except:
        # try trust-constr minimization first
        method = 'trust-constr'
        result = minimize(box_energy, x, (simulation.context, n), method='trust-constr', constraints=cons, jac=jacobian,
                          options={'disp': 1, 'xtol': 1e-08, 'gtol': 1e-08, 'maxiter': 1})

    # save and print the result
    x_new = result.x
    print("result", x_new)

    # update simulation with minimized positions and box vectors
    frac_positions_arr = np.empty([n, 3])
    for i in range(int(n)):
        frac_positions_arr[i][:] = x_new[i * 3:i * 3 + 3]

    # convert result to box vector
    box_parameter_new = [x_new[n * 3], x_new[n * 3 + 1], x_new[n * 3 + 2], x_new[n * 3 + 3], x_new[n * 3 + 4],
                         x_new[n * 3 + 5]]
    box_vectors_new = para2matrix(box_parameter_new, radians=False, format="lower")
    a = np.array([box_vectors_new[0][0], 0, 0])
    b = np.array([box_vectors_new[1][0], box_vectors_new[1][1], 0])
    c = np.array([box_vectors_new[2][0], box_vectors_new[2][1], box_vectors_new[2][2]])
    simulation.context.setPeriodicBoxVectors(a, b, c)

    # convert fractional position to atom position
    positions_arr = np.matmul(box_vectors_new, frac_positions_arr.T).T
    simulation.context.setPositions(positions_arr)

    # Calculate state and energy
    min_state = simulation.context.getState(getEnergy=True, getPositions=True, getForces=True)
    min_potential = min_state.getPotentialEnergy()
    simulation.saveState('data/final_states/' + cod_id + '_final.xml')
    print('Final Energy = ' + str(min_potential))
    simulation.step(1)

    # Load mdtraj data
    initial = mdtraj.load_pdb('data/PDB_supercell/' + cod_id + '_supercell.pdb')
    final = mdtraj.load_pdb('data/minimized_PDB_supercell/' + cod_id + '.pdb')

    ### Calculate RMSD
    # mdtraj RMSD
    rmsd_MDTraj = mdtraj.rmsd(initial, final)
    print("MDTraj rmsd", rmsd_MDTraj[0])
    # Tobias RMSD
    rmsd_Tobias = get_rmsd('data/PDB_supercell/' + cod_id + '_supercell.pdb',
                           'data/minimized_PDB_supercell/' + cod_id + '.pdb')
    print("Tobias rmsd", rmsd_Tobias)
    rmsd_Tobias = get_rmsd_n('data/PDB_supercell/' + cod_id + '_supercell.pdb',
                        'data/minimized_PDB_supercell/' + cod_id + '.pdb', "o")
    print("Tobias rmsd", rmsd_Tobias)

    # Remove Periodicity RMSD
    rmsd_r = get_rmsd_r('data/PDB_supercell/' + cod_id + '_supercell.pdb',
                        'data/minimized_PDB_supercell/' + cod_id + '.pdb')
    print("rmsd_r", rmsd_r)
    rmsd_r = get_rmsd_n('data/PDB_supercell/' + cod_id + '_supercell.pdb',
                        'data/minimized_PDB_supercell/' + cod_id + '.pdb', "r")
    print("rmsd_r", rmsd_r)
    # RMSD 20
    rmsd20 = get_rmsd_n('data/PDB_supercell/' + cod_id + '_supercell.pdb',
                        'data/minimized_PDB_supercell/' + cod_id + '.pdb', 20)
    print("rmsd20", rmsd20)

    data.append(
        {
            'COD ID': cod_id,
            'RMSD_mdtraj': rmsd_MDTraj[0],
            'RMSD_Tobias': rmsd_Tobias,
            'RMSD (Remove Periodicity)': rmsd_r,
            'RMSD20': rmsd20,
            'RMSD20 MAX': max(rmsd20),
            'RMSD20 MIN': min(rmsd20),
            'RMSD20 Mean': sum(rmsd20) / len(rmsd20),
            'method': method,
            'Original Energy': orig_potential,
            'Minimized Energy (OpenMM)': mid_potential,
            'Minimized Energy (Cell Minimization)': min_potential,
            'Original Energy/molecule': orig_potential / number_of_molecule_in_supercell,
            'Minimized Energy (OpenMM)/molecule': mid_potential / number_of_molecule_in_supercell,
            'Minimized Energy (Cell Minimization)/molecule': min_potential / number_of_molecule_in_supercell,
            'Original Box Vectors': A,
            'Minimized Box Vectors': np.array([a, b, c]),
            'Original Box parameters': box_parameter,
            'Minimized Box parameter': box_parameter_new,

        })
    with open('data/rmsd_values.txt', 'a') as f:
        f.write('%s\t%s\n' % (cod_id, sum(rmsd20) / len(rmsd20)))

# except Exception as e:
#    logging.error('Generic error with ID %s' % cod_id)
#    logging.error(e)
#    continue

try:
    d = pd.DataFrame(data)
    d.to_csv('data/minimization_results.csv')
    d.to_pickle('data/minimization_results.pkl')
except Exception as e:
    logging.error('Error with save of results data')
