# Import all COD files and generate all PDB files

import urllib.request
import os
from openbabel import pybel
import pymol


file = open('COD_ID_List.txt')

for line in file:
    url = 'http://www.crystallography.net/cod/' + line.strip() + '.cif'
    filename = 'data/CIF/' + line.strip() + '.cif'
    urllib.request.urlretrieve(url, filename)

for f in os.listdir('data/CIF'):
    mol = next(pybel.readfile("cif", 'data/CIF/' + f))
    id = f.split('.')[0]
    mol.write('pdb', 'data/PDB/' + id + '.pdb', overwrite=True)

