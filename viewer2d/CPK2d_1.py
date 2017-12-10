from Bio.PDB import *

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from molecular_data import color

#Desglosamos archivo PDB
parser = PDBParser(QUIET=True,PERMISSIVE=True)
structure = parser.get_structure('MACROH2A', 'data/1yd9.pdb')

#Creamos el modelo
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for atom in structure.get_atoms():
	x, y, z = atom.coord
	id = atom.get_id()
	ax.scatter(x, y, z, c=color(id), marker='o')
	ax.axis("off")

plt.show()