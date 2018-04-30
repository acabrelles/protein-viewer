from vispy import app,visuals,scene
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

Plot3D = scene.visuals.create_visual_node(visuals.LinePlotVisual)

canvas = scene.SceneCanvas(keys='interactive', app='pyqt4', bgcolor='white', title='Spheres',
                           size=(W,H), show=True)

view = canvas.central_widget.add_view()
view.camera = 'turntable'

Plot3D(coordinates, width=0.0, color='red', edge_color='b', symbol='o', face_color=(0.2, 0.2, 1, 1.0),parent=view.scene)

canvas.app.run()