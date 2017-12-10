from Bio.PDB import *

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from molecular_data import resdict, restype, color, colorrgba

#Desglosamos archivo PDB
parser = PDBParser(QUIET=True,PERMISSIVE=True)
structure = parser.get_structure('MACROH2A', 'data/1yd9.pdb')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#Seleccionamos los residuos, solo contando los aminoacidos
for resname, residuetype in resdict.iteritems():
    residues = [residue for residue in structure.get_residues() if residue.get_resname() == resname]
    rescoord = []
    color = colorrgba(residuetype)
    
    for residue in residues:
        atoms = [atom for atom in residue.get_atoms()]
        coordinates = [atom.coord for atom in atoms]
        rescoord.append(np.array(coordinates))
        
    if len(rescoord)>1:
        rescoord = np.concatenate(rescoord)
        
    if not len(residues)==0:
        x, y, z =zip(*rescoord)
        ax.scatter(x, y, z, c=color, marker='o')
        ax.axis("off")

#Seleccionamos el resto de 'residuos' que ha determinado el parser, sin incluir aguas
residues_1 = [residue for residue in structure.get_residues() if residue.get_resname() != 'HOH']
residues_2 = [residue for residue in structure.get_residues() if residue.get_resname() in resdict.keys()]

residues_pink = list(set(residues_1)-set(residues_2))
rescoordpink = []

for residue in residues_pink:
    atomspink = [atom for atom in residue.get_atoms()]
    coordinatespink = [atom.coord for atom in atomspink]
    rescoordpink.append(np.array(coordinatespink))
    
if len(rescoordpink)>1:
    rescoordpink = np.concatenate(rescoordpink)
    
xp, yp, zp = zip(*rescoordpink)
ax.scatter(xp, yp, zp, c='pink', marker='o')
ax.axis("off")

plt.show()