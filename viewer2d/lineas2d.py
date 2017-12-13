from Bio.PDB import *

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from molecular_data import colors

#Desglosamos archivo PDB
parser = PDBParser(QUIET=True,PERMISSIVE=True)
structure = parser.get_structure('MACROH2A', 'data/1yd9.pdb')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for chain in structure.get_chains():
        can_atoms = [atom for atom in chain.get_atoms() if atom.get_name() == 'CA' or atom.get_name() == 'N']
        can_coordinates = [atom.coord for atom in can_atoms]
        x,y,z=zip(*can_coordinates)
        ccolor = np.random.rand(3,1)
        ax.plot(x, y, z, c=ccolor, linewidth=2)
        ax.scatter(x, y, z, c=ccolor, marker='o')
        ax.axis("off")
        
plt.show()
