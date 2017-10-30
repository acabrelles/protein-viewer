from Bio.PDB import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

parser = PDBParser()
structure = parser.get_structure('MACROH2A', '1yd9.pdb')

chain_list=Selection.unfold_entities(structure,'C')
residue_list=Selection.unfold_entities(chain_list[0], 'R')
atom_list=Selection.unfold_entities(residue_list,'A')

coord=[]
for atom in atom_list:
    if atom.id == 'CA' or atom.id == 'N': 
        coord.append(atom.get_vector())
x,y,z=zip(*coord)

fig = plt.figure()
ax=Axes3D(fig)
ax.plot(x,y,z,linewidth=3)
ax.axis("off")

plt.show()