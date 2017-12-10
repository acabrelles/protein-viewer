from Bio.PDB import *

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from molecular_data import colors

#Desglosamos archivo PDB
parser = PDBParser(QUIET=True,PERMISSIVE=True)
structure = parser.get_structure('MACROH2A', 'data/1yd9.pdb')

#Creamos el modelo
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for atom_type, color in colors.iteritems():
        atoms = [atom for atom in structure.get_atoms() if atom.get_id() == atom_type]
        coordinates = [atom.coord for atom in atoms]
        if not len(atoms)==0:
            x,y,z=zip(*coordinates)
            ax.scatter(x, y, z, c=color, marker='o')
            ax.axis("off")


#Creamos un solo array para todos los atomos que no se han representado ya

atoms_1 = [atom for atom in structure.get_atoms()]
atoms_2 = [atom for atom in structure.get_atoms() if atom.get_id() in colors.keys()]

atoms_pink = list(set(atoms_1)-set(atoms_2))
coordinates_pink = [atom.coord for atom in atoms_pink]
xp, yp, zp = zip(*coordinates_pink)

ax.scatter(xp, yp, zp, c='pink', marker='o')
ax.axis("off")

plt.show()