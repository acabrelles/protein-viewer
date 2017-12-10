from Bio.PDB import *

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from molecular_data import resdict, restype, color, colorrgba

#Desglosamos archivo PDB
parser = PDBParser(QUIET=True,PERMISSIVE=True)
structure = parser.get_structure('MACROH2A', 'data/1yd9.pdb')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#Seleccionamos los residuos, solo contando los aminoacidos
residues = [residue for residue in structure.get_residues() if residue.get_resname() in resdict.keys()]
for residue in residues:
    resid = residue.get_resname()
    rescolor = colorrgba(restype(resid))
    atoms = [atom for atom in residue.get_atoms()]
    for atom in atoms:
        x, y, z = atom.coord
        ax.scatter(x, y, z, c=rescolor, marker='o')
        ax.axis("off")

#Seleccionamos el resto de 'residuos' que ha determinado el parser, sin incluir aguas
residues2 = [residue for residue in structure.get_residues() if not residue in residues and residue.get_resname() != 'HOH']
for residue in residues2:
    atoms = [atom for atom in residue.get_atoms()]
    for atom in atoms:
        x, y, z = atom.coord
        ax.scatter(x, y, z, c='pink', marker='o')
        ax.axis("off")

plt.show()