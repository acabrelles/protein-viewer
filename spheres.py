from vispy import scene
from vispy.visuals.transforms import STTransform

import numpy as np

from Bio.PDB import PDBParser,DSSP

from molecular_data import crgbaDSSP, restype, colorrgba, vrad, resdict

pdbdata = 'data/1yd9.pdb'
parser = PDBParser(QUIET=True, PERMISSIVE=True)
structure = parser.get_structure('model',pdbdata)

def centroid(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return sum_x/length, sum_y/length, sum_z/length

atoms = [atom for atom in structure.get_atoms()]
natoms = len(atoms)
#atom coordinates
coordinates = np.array([atom.coord for atom in atoms])
center = centroid(coordinates)
coordinates -= center
#atom color
color = [colorrgba(atom.get_id()) for atom in atoms]
#atom radius
radius = [vrad(atom.get_id()) for atom in atoms]

W,H = 1200, 800

canvas = scene.SceneCanvas(keys='interactive', app='pyqt4', bgcolor='white',
                           size=(800, 600), show=True)

view = canvas.central_widget.add_view()
view.camera = 'arcball'

sphere_list = []

for i in range(len(atoms)):
    sphere = scene.visuals.Sphere(radius=radius[i], method='ico', parent = view.scene, color=color[i])
    sphere_list.append(sphere)
for i in range(len(coordinates)):
    sphere_list[i].transform = STTransform(translate=coordinates[i])