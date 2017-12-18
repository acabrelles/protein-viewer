from Bio.PDB import PDBParser,DSSP

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from molecular_data import colors, colorsDSSP, resdict, restype, colorrgba

import sys, os

if len(sys.argv[1:])!=2:
        raise SystemExit("Uso: %s Archivo_PDB Visualizacion(cpk, backbone, aminoacid, dssp)" % os.path.basename(sys.argv[0]))

#Desglosamos archivo PDB
pdbdata = sys.argv[1]
parser = PDBParser(QUIET=True,PERMISSIVE=True)
structure = parser.get_structure('model', pdbdata)
        
#Hacemos la prediccion DSSP
model = structure[0]
dssp = DSSP(model, pdbdata)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

def cpk2d():
    for atom_type, color in colors.iteritems():
            atoms = [atom for atom in structure.get_atoms() if atom.get_id() == atom_type]
            coordinates = [atom.coord for atom in atoms]
            if not len(atoms)==0:
                x, y, z=zip(*coordinates)
                ax.scatter(x, y, z, c=color, marker='o')
                
    #Seleccionamos los atomos no identificados
    atoms_1 = [atom for atom in structure.get_atoms()]
    atoms_2 = [atom for atom in structure.get_atoms() if atom.get_id() in colors.keys()]
    atoms_pink = list(set(atoms_1)-set(atoms_2))
    
    coordinates_pink = [atom.coord for atom in atoms_pink]
    xp, yp, zp = zip(*coordinates_pink)
    ax.scatter(xp, yp, zp, c='pink', marker='o')
    
    ax.axis("off")
    plt.show()
    
def bb2d():

    for chain in structure.get_chains():
        can_atoms = [atom for atom in chain.get_atoms() if atom.get_name() == 'CA' or atom.get_name() == 'N']
        can_coordinates = [atom.coord for atom in can_atoms]
        x,y,z=zip(*can_coordinates)
        ccolor = np.random.rand(3,1)
        ax.plot(x, y, z, c=ccolor, linewidth=2)
        ax.scatter(x, y, z, c=ccolor, marker='o')
        
    ax.axis("off")    
    plt.show()
    
def aa2d():
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
    
def dssp2d():
    ax = fig.add_subplot(111, projection='3d')
    
    #Creamos las listas de residuos vinculadas a su prediccion
    residues = [residue for residue in structure.get_residues() if residue.get_resname() in resdict.keys()]
    struct3 = [dssp[key][2] for key in list(dssp.keys())]
    respred = zip(struct3,residues)

    #Creamos la nube de puntos por prediccion
    for prediction, color in colorsDSSP.iteritems():
        residuesp = [residue[1] for residue in respred if residue[0] == prediction]
        predcoord_can = []
        for residue in residuesp:
            atomsp = [atom for atom in residue.get_atoms() if atom.get_name() == 'CA' or atom.get_name() == 'N']
            coordinatesp = [atom.coord for atom in atomsp]
            predcoord_can.append(np.array(coordinatesp))
        if len(predcoord_can)>1:
            predcoord_can = np.concatenate(predcoord_can)
        if not len(residuesp)==0:
            x, y, z = zip(*predcoord_can)
            ax.scatter(x, y, z, c=color, marker='o')

    #Creamos las cadenas que unen los atomos
    for chain in structure.get_chains():
        can_atoms = [atom for atom in chain.get_atoms() if atom.get_name() == 'CA' or atom.get_name() == 'N']
        can_coordinates = [atom.coord for atom in can_atoms]
        x, y, z = zip(*can_coordinates)
        ccolor = np.random.rand(3,1)
        ax.plot(x, y, z, c=ccolor, linewidth=1)
    
    ax.axis("off")
    plt.show()

if sys.argv[2]=='cpk':
    cpk2d()
elif sys.argv[2]=='backbone':
    bb2d()
elif sys.argv[2]=='aminoacid':
    aa2d()
elif sys.argv[2]=='dssp':
    dssp2d()
else:
    print 'Not recognized visualization mode %s' % sys.argv[2]