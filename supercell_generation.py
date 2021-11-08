# Generate 3x3x3 supercells for all pdb files and apply cryst1 modifications

import pymol
import os
from openbabel import pybel

for f in os.listdir('data/PDB'):
    print(f)
    # reinitialize to clear pymol
    pymol.cmd.reinitialize()
    pymol.cmd.load('data/PDB/' + f)
    id = f.split('.')[0]
    # load supercell and run for molecule
    pymol.cmd.do('run https://raw.githubusercontent.com/speleo3/pymol-psico/master/psico/xtal.py')
    pymol.cmd.do('supercell 3,3,3,' + id)
    pymol.cmd.delete(id)
    # Make the supercell file if needed
    path = '/home/qualenal/Scripts/OpenFFCrystalBenchmarking/data/PDB_supercell/'+id+'_supercell.pdb'
    try:
        os.mknod(path)
    except FileExistsError:
        pass
    pymol.cmd.save(path,format='pdb',quiet='0')
    try:
        mol = next(pybel.readfile('pdb',path))
        mol.unitcell.SetData(mol.unitcell.GetA()*3, mol.unitcell.GetB()*3, mol.unitcell.GetC()*3, mol.unitcell.GetAlpha(), mol.unitcell.GetBeta(), mol.unitcell.GetGamma())
        mol.write('pdb',path,overwrite=True)
    except AttributeError:
        print(id)


