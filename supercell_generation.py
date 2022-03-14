# Generate supercells for all pdb files and apply cryst1 symmetry operations
import math
import pymol
import os
from openbabel import pybel

for f in os.listdir('data/PDB'):
    print(f)
    # load pdb file to get unit cell dimensions and pick # of cells to ensure each dimension is >18 A
    try:
        mol = next(pybel.readfile('pdb','data/PDB/' + f))
    except Exception as e:
        print(e)
        continue
    # Extract cell length
    a_mult = math.ceil(24.0/mol.unitcell.GetA())
    b_mult = math.ceil(24.0/mol.unitcell.GetB())
    c_mult = math.ceil(24.0/mol.unitcell.GetC())
    # add one extra cell in the smallest direction
    if a_mult < b_mult and a_mult < c_mult:
        a_mult = a_mult + 1
    elif b_mult < c_mult and b_mult < a_mult:
        b_mult = b_mult + 1
    elif c_mult < a_mult and c_mult < b_mult:
        c_mult = c_mult + 1
    # reinitialize to clear pymol
    pymol.cmd.reinitialize()
    pymol.cmd.load('data/PDB/' + f)
    id = f.split('.')[0]
    # load supercell and run for molecule
    pymol.cmd.do('run https://raw.githubusercontent.com/speleo3/pymol-psico/master/psico/xtal.py')
    pymol.cmd.do('supercell %s,%s,%s,%s' % (a_mult, b_mult, c_mult, id))
    pymol.cmd.delete(id)
    # Make the supercell file if needed
    path = './data/PDB_supercell/'+id+'_supercell.pdb'
    try:
        os.mknod(path)
    except FileExistsError:
        pass
    pymol.cmd.save(path,format='pdb',quiet='0')
    try:
        mol_sc = next(pybel.readfile('pdb', path))
        mol_sc.unitcell.SetData(mol_sc.unitcell.GetA() * a_mult, mol_sc.unitcell.GetB() * b_mult, mol_sc.unitcell.GetC() * c_mult, mol_sc.unitcell.GetAlpha(), mol_sc.unitcell.GetBeta(), mol_sc.unitcell.GetGamma())
        mol_sc.write('pdb', path, overwrite=True)
    except AttributeError:
        print(id)


