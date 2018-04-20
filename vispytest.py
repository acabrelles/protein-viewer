from Bio.PDB import PDBParser,DSSP

from vispy import gloo, app, visuals

from vispy.util.transforms import perspective, translate

import numpy as np

from molecular_data import crgbaDSSP, restype, colorrgba, vrad, resdict

W,H = 1200, 800

vertex= """
#version 120

uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_projection;

attribute vec3  a_position;
attribute vec4  a_color;
attribute float a_radius;

varying vec4 v_color;
varying float v_radius;
varying v_eye_position;

void main (void) {
    v_radius = a_radius;
    v_color = a_color;
    
    v_eye_position = u_view * u_model * vec4(a_position,1.0);
    
    gl_Position = u_projection * v_eye_position;
    gl_PointSize = v_radius
}

"""
fragment= """
varying vec4 v_color;

void main (void) {
    gl_FragColor = v_color;
}
"""

#Function to define the center of the point cloud
def centroid(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return sum_x/length, sum_y/length, sum_z/length

class Canvas(app.Canvas):
    
    visualization_modes = ['cpk','backbone','aminoacid','dssp']
    
    def __init__(self, pdbdata, mode='cpk'):
        
        #Startup
        app.Canvas.__init__(self, keys='interactive', size =(W,H))
        
        #Loading shaders
        self.program = gloo.Program(vertex, fragment)
        
        #Analyze pdb file
        self.parser = PDBParser(QUIET=True, PERMISSIVE=True)
        self.structure = self.parser.get_structure('model',pdbdata)
        
        #DSSP prediction
        self.model = self.structure[0]
        self.dssp = DSSP(self.model, pdbdata)
        
        #Mode selection
        if mode not in Canvas.visualization_modes:
            raise Exception('Not recognized visualization mode %s' % mode)
        self.mode = mode
        
        #Camera settings
        self.translate = 120
        self.view = translate((0,0, -self.translate), dtype=np.float32)
        self.model = np.eye(4, dtype=np.float32)
        self.projection = perspective(45.0, self.size[0] / float(self.size[1]), 1.0, 1000.0)
        
        self.program['u_model'] = self.model
        self.program['u_view'] = self.view
        self.program['u_projection'] = self.projection
        
        #Load data depending on the mdoe
        self.atom_information()
        self.load_data()
        
        
        
    def atom_information(self):
        
        """Determines the coordinates, colors and sizes of the atoms depending on the mode"""
        
        if self.mode == 'cpk':
            #list of atoms
            self.atoms = [atom for atom in self.structure.get_atoms()]
            self.natoms = len(self.atoms)
            #atom coordinates
            self.coordinates = np.array([atom.coord for atom in self.atoms])
            self.center = centroid(self.coordinates)
            self.coordinates -= self.center
            #atom color
            self.color = [colorrgba(atom.get_id()) for atom in self.atoms]
            #atom radius
            self.radius = [vrad(atom.get_id() for atom in self.atoms)]
            
        elif self.mode == 'aminoacid':
            #list of atoms
            self.atoms = [atom for atom in self.structure.get_atoms()]
            self.natoms = len(self.atoms)
            #atom coordinates
            self.coordinates = np.array([atom.coord for atom in self.atoms])
            self.center = centroid(self.coordinates)
            self.coordinates -= self.center
            #atom color
            self.color = [colorrgba(restype(atom.get_parent().resname)) for atom in self.atoms]
            #atom radius
            self.radius = [vrad(atom.get_id() for atom in self.atoms)]
            
        elif self.mode == 'backbone':
            #list of atoms
            self.atoms = [atom for atom in self.structure.get_atoms() if atom.get_name() == 'CA' or atom.get_name() == 'N']
            self.natoms = len(self.atoms)
            #atom coordinates
            self.coordinates = np.array([atom.coord for atom in self.atoms])
            self.center = centroid(self.coordinates)
            self.coordinates -= self.center
            #atom color
            self.color = []
            self.chains = []
            for chain in self.structure.get_chains():
                self.chains.append(chain)
                self.chain_length = len([atom for atom in chain.get_atoms() if atom.get_name() =='CA' or atom.get_name() =='N'])
                self.chain_color = np.append(np.random.rand(1,3),[1.0])
                self.color.append(np.tile(self.chain_color,(self.chain_length,1)))
            if len(self.chains)>1:
                self.color = np.concatenate(self.color)
            #atom radius
            self.radius = [vrad(atom.get_id() for atom in self.atoms)]
                
        elif self.mode == 'dssp':
            #list of atoms
            self.atoms = [atom for atom in self.structure.get_atoms() if atom.get_name() == 'CA' or atom.get_name() == 'N']
            self.natoms = len(self.atoms)
            #atom coordinates
            self.coordinates = np.array([atom.coord for atom in self.atoms])
            self.center = centroid(self.coordinates)
            self.coordinates -= self.center
            #atom color
            self.struct3 = [self.dssp[key][2] for key in list(self.dssp.keys())]
            self.residues = [residue for residue in self.structure.get_residues() if residue.get_resname() in resdict.keys()]
            self.color = []
            for i in range(len(self.struct3)):
                self.dsspcolor = crgbaDSSP(self.struct3[i])
                self.natoms = len([atom for atom in self.residues[i] if atom.get_name() =='CA' or atom.get_name() == 'N'])
                self.color.append(np.tile(self.dsspcolor,(self.natoms,1)))
            if len(self.struct3)>1:
                self.color = np.concatenate(self.color)
            #atom radius
            self.radius = [vrad(atom.get_id() for atom in self.atoms)]
            
    def load_data(self):
        
        """Make an array with all the data and load it into VisPy Gloo"""
        
        data = np.zeros(self.natoms, [('a_position', np.float32, 3),
                            ('a_color', np.float32, 4),
                            ('a_radius', np.float32, 1)])

        data['a_position'] = self.coordinates
        data['a_color'] = self.color
        data['a_radius'] = self.radius*self.pixel_scale

        self.program.bind(gloo.VertexBuffer(data))
    
    def on_resize(self, event):
        gloo.set_viewport(0, 0, event.physical_size[0], event.physical_size[1])
        self.projection = perspective(45.0, event.size[0] / float(event.size[1]), 1.0, 1000.0)
        self.program['u_projection'] = self.projection
    
    def on_draw(self,event):
        gloo.clear()
        self.program.draw('points')

pdbdata = 'data/1yd9.pdb'
mvc = Canvas(pdbdata, mode='cpk')
app.run()