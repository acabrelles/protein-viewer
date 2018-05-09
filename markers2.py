from vispy import app,visuals,scene
from vispy.visuals.transforms import STTransform
from vispy.scene import ArcballCamera

import numpy as np

from Bio.PDB import PDBParser,DSSP

from molecular_data import crgbaDSSP, restype, colorrgba, vrad, resdict

pdbdata = 'data/1yd9.pdb'
parser = PDBParser(QUIET=True, PERMISSIVE=True)
structure = parser.get_structure('model',pdbdata)

pmodel = structure[0]
dssp = DSSP(pmodel, pdbdata)

def centroid(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return sum_x/length, sum_y/length, sum_z/length

#list of atoms
atoms = [atom for atom in structure.get_atoms() if atom.get_name() == 'CA' or atom.get_name() == 'N']
natoms = len(atoms)
#atom coordinates
coordinates = np.array([atom.coord for atom in atoms])
center = centroid(coordinates)
coordinates -= center
#atom color
color = []
chains = []
chain_coords = []
chain_colors = []
chain_radius = []
#chains
for chain in structure.get_chains():
    chains.append(chain)
    can_coord = np.array([atom.coord for atom in chain.get_atoms() if atom.get_name() =='CA' or atom.get_name() =='N'])
    can_coord -= center
    chain_coords.append(can_coord)
    chain_length = len(can_coord)
    chain_color = np.random.rand(1,3)
    chain_colors.append(chain_color)
    color.append(np.tile(chain_color,(chain_length,1)))
    chain_radius.append([vrad(atom.get_id()) for atom in chain.get_atoms() if atom.get_name() =='CA' or atom.get_name() =='N']) 
if len(chains)>1:
    color = np.concatenate(color)
#dssp
color2=[]
struct3 = [dssp[key][2] for key in list(dssp.keys())]
residues = [residue for residue in structure.get_residues() if residue.get_resname() in resdict.keys()]
for i in range(len(struct3)):
    dsspcolor = crgbaDSSP(struct3[i])[0:3]
    n_atoms = len([atom for atom in residues[i] if atom.get_name() =='CA' or atom.get_name() == 'N'])
    color2.append(np.tile(dsspcolor,(n_atoms,1)))
if len(struct3)>1:
    color2 = np.concatenate(color2)
#atom radius
radius = [4*vrad(atom.get_id()) for atom in atoms if atom.get_name() =='CA' or atom.get_name() == 'N']

W,H = 1200, 800

Plot3D = scene.visuals.create_visual_node(visuals.LinePlotVisual)

canvas = scene.SceneCanvas(keys='interactive', app='pyqt4', bgcolor='white', title='Spheres',
                           size=(W,H), show=True)

view = canvas.central_widget.add_view()
view.camera = ArcballCamera(fov=95, distance=(max(abs(np.concatenate(coordinates))) + 40))

print color2

Plot3D(coordinates, marker_size=radius, width=0.0, color='red', face_color=color2, edge_color='b', symbol='o', parent=view.scene)
for i in range(len(chains)):
    Plot3D(chain_coords[i], marker_size=0.0, width=10.0, color=chain_colors[i], parent=view.scene)

canvas.app.run()