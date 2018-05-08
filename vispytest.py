from Bio.PDB import PDBParser

from vispy import gloo, app, visuals

from vispy.util.transforms import perspective, translate, rotate


import numpy as np

from molecular_data import crgbaDSSP, restype, colorrgba, vrad, resdict
from vispy.util.quaternion import Quaternion


W,H = 1200, 800

vertex = """
#version 120

uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_projection;
uniform vec3 u_light_position;
uniform vec3 u_light_spec_position;

attribute vec3  a_position;
attribute vec3  a_color;
attribute float a_radius;

varying vec3  v_color;
varying vec4  v_eye_position;
varying float v_radius;
varying vec3  v_light_direction;

void main (void) {
    v_radius = a_radius;
    v_color = a_color;

    v_eye_position = u_view * u_model * vec4(a_position,1.0);
    v_light_direction = normalize(u_light_position);
    float dist = length(v_eye_position.xyz);

    gl_Position = u_projection * v_eye_position;

    vec4  proj_corner = u_projection * vec4(a_radius, a_radius, v_eye_position.z, v_eye_position.w);  // # noqa
    gl_PointSize = 512.0 * proj_corner.x / proj_corner.w;
}
"""

fragment = """
#version 120

uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_projection;
uniform vec3 u_light_position;
uniform vec3 u_light_spec_position;

varying vec3  v_color;
varying vec4  v_eye_position;
varying float v_radius;
varying vec3  v_light_direction;
void main()
{
    vec2 texcoord = gl_PointCoord* 2.0 - vec2(1.0);
    float x = texcoord.x;
    float y = texcoord.y;
    float d = 1.0 - x*x - y*y;
    if (d <= 0.0)
        discard;

    float z = sqrt(d);
    vec4 pos = v_eye_position;
    pos.z += v_radius*z;
    vec3 pos2 = pos.xyz;
    pos = u_projection * pos;
    // gl_FragDepth = 0.5*(pos.z / pos.w)+0.5;
    vec3 normal = vec3(x,y,z);
    float diffuse = clamp(dot(normal, v_light_direction), 0.0, 1.0);

    // Specular lighting.
    vec3 M = pos2.xyz;
    vec3 O = v_eye_position.xyz;
    vec3 L = u_light_spec_position;
    vec3 K = normalize(normalize(L - M) + normalize(O - M));
    // WARNING: abs() is necessary, otherwise weird bugs may appear with some
    // GPU drivers...
    float specular = clamp(pow(abs(dot(normal, K)), 40.), 0.0, 1.0);
    vec3 v_light = vec3(1., 1., 1.);
    gl_FragColor.rgb = (.15*v_color + .55*diffuse * v_color
                        + .35*specular * v_light);
}
"""

#Function to define the center of the point cloud
def centroid(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return sum_x/length, sum_y/length, sum_z/length

#Function to convert x,y coordinates to Quaternion for camera
def _arcball(x, y, w, h):
        """Convert x,y coordinates to w,x,y,z Quaternion parameters
        Adapted from:
        linalg library
        Copyright (c) 2010-2015, Renaud Blanch <rndblnch at gmail dot com>
        Licence at your convenience:
        GPLv3 or higher <http://www.gnu.org/licenses/gpl.html>
        BSD new <http://opensource.org/licenses/BSD-3-Clause>
        """
        x = float(x)
        y = float(y)
        w = float(w)
        h = float(h)
        r = (w + h) / 2.
        x, y = -(2. * x - w) / r, -(2. * y - h) / r
        h = np.sqrt(x*x + y*y)
        return (0., -x/h, y/h, 0.) if h > 1. else (0., -x, y, np.sqrt(1. - h*h))

class Canvas(app.Canvas):
    
    visualization_modes = ['cpk','backbone','aminoacid', 'dssp']
    
    def __init__(self, pdbdata, mode='cpk'):
        
        #Startup
        app.Canvas.__init__(self, keys='interactive', size =(W,H))
        
        #Loading shaders
        self.program = gloo.Program(vertex, fragment)
        
        #Analyze pdb file
        self.parser = PDBParser(QUIET=True, PERMISSIVE=True)
        self.structure = self.parser.get_structure('model',pdbdata)

        #Mode selection
        if mode not in Canvas.visualization_modes:
            raise Exception('Not recognized visualization mode %s' % mode)
        self.mode = mode
        
        #Camera settings
        self.translate = 100
        self.translate = max(-1, self.translate)
        self.view = translate((0,0, -self.translate), dtype=np.float32)

        self.model = np.eye(4, dtype=np.float32)
        self.projection = np.eye(4, dtype=np.float32)
        self.program['u_projection'] = self.projection
        
        self.quaternion = Quaternion()

        #Load data depending on the mdoe

        self.apply_zoom()
        self.atom_information()
        self.load_data()

        gloo.set_state(depth_test=True, clear_color='black')

        self.show()

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
            self.color = [np.array(colorrgba(atom.get_id())[0:3]) for atom in self.atoms]
            #atom radius
            self.radius = [vrad(atom.get_id()) for atom in self.atoms]
            
        elif self.mode == 'aminoacid':
            #list of atoms
            self.atoms = [atom for atom in self.structure.get_atoms() if atom.get_parent().resname != 'HOH']
            self.natoms = len(self.atoms)
            #atom coordinates
            self.coordinates = np.array([atom.coord for atom in self.atoms])
            self.center = centroid(self.coordinates)
            self.coordinates -= self.center
            #atom color
            self.color = [colorrgba(restype(atom.get_parent().resname))[0:3] for atom in self.atoms]
            #atom radius
            self.radius = [vrad(atom.get_id()) for atom in self.atoms]
            
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
                self.chain_color = np.random.rand(1,3)
                self.color.append(np.tile(self.chain_color,(self.chain_length,1)))
            if len(self.chains)>1:
                self.color = np.concatenate(self.color)
            #atom radius
            self.radius = [vrad(atom.get_id()) for atom in self.atoms]

    def load_data(self):
        
        """Make an array with all the data and load it into VisPy Gloo"""
        
        data = np.zeros(self.natoms, [('a_position', np.float32, 3),
                            ('a_color', np.float32, 3),
                            ('a_radius', np.float32, 1)])

        data['a_position'] = self.coordinates
        data['a_color'] = self.color
        data['a_radius'] = self.radius #*self.pixel_scale

        self.program.bind(gloo.VertexBuffer(data))

        self.program['u_model'] = self.model
        self.program['u_view'] = self.view
        self.program['u_light_position'] = 0., 0., 2.
        self.program['u_light_spec_position'] = -5., 5., -5.

        print 'Data loaded'
    
    def on_resize(self, event):
        width, height = event.size

    def apply_zoom(self):
        width, height = self.physical_size
        gloo.set_viewport(0, 0, width, height)
        self.projection = perspective(95.0, width / float(height), 1.0, 400.0)
        self.program['u_projection'] = self.projection
    
    def on_draw(self,event):
        gloo.clear()
        self.program.draw('points')

    def on_mouse_move(self,event):
        if event.button == 1 and event.last_event is not None:
            x0, y0 = event.last_event.pos
            x1, y1 = event.pos
            w, h = self.size
            self.quaternion = (self.quaternion *
                               Quaternion(*_arcball(x0, y0, w, h)) *
                               Quaternion(*_arcball(x1, y1, w, h)))
            self.model = self.quaternion.get_matrix()
            self.program['u_model'] = self.model
            self.update()
        elif event.button == 2 and event.last_event is not None:
            x0, y0 = event.last_event.pos
            x1, y1 = event.pos
            self.translate += (y1-y0)
            self.translate = max(-1, self.translate)
            self.view = translate((0,0, -self.translate), dtype=np.float32)
            self.program['u_view'] = self.view
            self.update()


pdbdata = 'data/1yd9.pdb'
mvc = Canvas(pdbdata, mode='aminoacid')
app.run()
