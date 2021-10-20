from openbabel import pybel
import pymol

mol = next(pybel.readfile("cif", 'data/CIF/1100249.cif'))

mol.write('pdb', 'data/PDB/1100249.pdb')

#pymol.cmd.do('run https://raw.githubusercontent.com/speleo3/pymol-psico/master/psico/xtal.py')
