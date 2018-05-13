#VisPy Spheres Class

import numpy as np
from vispy import app, gloo, visuals, scene

vertex = """
#version 120

uniform vec3 u_light_position;
uniform vec3 u_light_spec_position;

attribute vec3  a_position;
attribute vec3  a_color;
attribute float a_radius;

varying vec3  v_color;
varying vec4  v_eye_direction;
varying float v_radius;
varying vec3  v_light_direction;

varying float v_depth;
varying float v_depth_radius;

void main (void) {
    vec4 atom_pos = vec4(a_position, 1);
    
    // First decide where to draw this atom on screen
    vec4 fb_pos = $visual_to_framebuffer(atom_pos);
    gl_Position = $framebuffer_to_render(fb_pos);
    
    // Measure the orientation of the framebuffer coordinate system relative
    // to the atom
    vec4 x = $framebuffer_to_visual(fb_pos + vec4(100, 0, 0, 0));
    x = (x/x.w - atom_pos) / 100;
    vec4 z = $framebuffer_to_visual(fb_pos + vec4(0, 0, -100, 0));
    z = (z/z.w - atom_pos) / 100;
    
    // Use the x axis to measure radius in framebuffer pixels
    // (gl_PointSize uses the framebuffer coordinate system)
    vec4 radius = $visual_to_framebuffer(atom_pos + normalize(x) * a_radius);
    radius = radius/radius.w - fb_pos/fb_pos.w;
    gl_PointSize = length(radius);
    
    // Use the z axis to measure position and radius in the depth buffer
    v_depth = gl_Position.z / gl_Position.w;
    // gl_FragDepth uses the "render" coordinate system.
    vec4 depth_z = $framebuffer_to_render($visual_to_framebuffer(atom_pos + normalize(z) * a_radius));
    v_depth_radius = v_depth - depth_z.z / depth_z.w;
    
    v_light_direction = normalize(u_light_position);
    v_radius = a_radius;
    v_color = a_color;
}
"""

fragment = """
#version 120

varying vec3  v_color;
varying float v_radius;
varying vec3  v_light_direction;
varying float v_depth;
varying float v_depth_radius;

void main()
{
    // calculate xyz position of this fragment relative to radius
    vec2 texcoord = gl_PointCoord * 2.0 - vec2(1.0);
    float x = texcoord.x;
    float y = texcoord.y;
    float d = 1.0 - x*x - y*y;
    if (d <= 0.0)
        discard;
    float z = sqrt(d);
    vec3 normal = vec3(x,y,z);
    
    // Diffuse color
    float ambient = 0.3;
    float diffuse = dot(v_light_direction, normal);
    // clamp, because 0 < theta < pi/2
    diffuse = clamp(diffuse, 0.0, 1.0);
    vec3 light_color = vec3(1, 1, 1);
    vec3 diffuse_color = ambient + light_color * diffuse;

    // Specular color
    //   reflect light wrt normal for the reflected ray, then
    //   find the angle made with the eye
    vec3 eye = vec3(0, 0, -1);
    float specular = dot(reflect(v_light_direction, normal), eye);
    specular = clamp(specular, 0.0, 1.0);
    // raise to the material's shininess, multiply with a
    // small factor for spread
    specular = pow(specular, 80);
    vec3 specular_color = light_color * specular;
    
    gl_FragColor = vec4(v_color * diffuse_color + specular_color, 1);
    gl_FragDepth = v_depth - .5 * z * v_depth_radius;
}
"""

class SpheresVisual(visuals.Visual):
    """Visual that draws many spheres.
    
    Parameters
    ----------
    coordinates: array of coordinates
    color: array of colors
    radius: array of radius
    """
    def __init__(self, coordinates, color, radius):
        visuals.Visual.__init__(self, vertex, fragment)
        
        self.natoms = len(coordinates)
        
        #Loading data and type
        self._load_data()
        self._draw_mode = 'points'
        self.set_gl_state('translucent', depth_test=True, cull_face=False)        
        
    def _load_data(self):
        """Make an array with all the data and load it into VisPy Gloo"""
        data = np.zeros(self.natoms, [('a_position', np.float32, 3),
                            ('a_color', np.float32, 4),
                            ('a_radius', np.float32, 1)])

        data['a_position'] = coordinates
        data['a_color'] = color
        data['a_radius'] = radius#*view.transforms.pixel_scale

        self.shared_program.bind(gloo.VertexBuffer(data))
        
        self.shared_program['u_light_position'] = 5., -5., 5.
    
    def _prepare_transforms(self,view):
        view.view_program.vert['visual_to_framebuffer'] = view.get_transform('visual', 'framebuffer')
        view.view_program.vert['framebuffer_to_visual'] = view.get_transform('framebuffer', 'visual')
        view.view_program.vert['framebuffer_to_render'] = view.get_transform('framebuffer', 'render')

#VisPyViewer

from vispy import scene
from Bio.PDB import PDBParser,DSSP
from molecular_data import crgbaDSSP, restype, colorrgba, vrad, resdict
        
def centroid(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return sum_x/length, sum_y/length, sum_z/length
        
class VisPyViewer(object):
    
    visualization_modes = ['cpk','backbone','aminoacid', 'dssp']
    
    def __init__(self, pdbdata, mode='cpk'):
        
        #Analyze pdb file
        self.parser = PDBParser(QUIET=True, PERMISSIVE=True)
        self.structure = self.parser.get_structure('model',pdbdata)
        
        #DSSP prediction
        self.pmodel = self.structure[0]
        self.dssp = DSSP(self.pmodel, pdbdata)
        
        #Mode selection
        if mode not in VisPyViewer.visualization_modes:
            raise Exception('Not recognized visualization mode %s' % mode)
        self.mode = mode
        
        #Data selection
        self.atom_information()
        print len(self.coordinates)
        print len(self.color)
        print len(self.radius)
        self.radius = max(abs(np.concatenate(self.coordinates)))
        
        #Canvas + camera
        canvas = scene.SceneCanvas(keys='interactive', app='pyqt4', bgcolor='white', size=(1200,800), show=True)
        view = canvas.central_widget.add_view()
        view.camera = scene.ArcballCamera(fov=70, distance=(self.radius+40))
        
        #Load visual and apply it
        Spheres = scene.visuals.create_visual_node(SpheresVisual)
        vis_atoms=[Spheres(coordinates,color,radius,parent=view.scene)]
        
        #Run the program
        canvas.app.run()
        
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
            self.radius = np.array([vrad(atom.get_id()) for atom in self.atoms])
            
        elif self.mode == 'aminoacid':
            #list of atoms
            self.atoms = [atom for atom in self.structure.get_atoms() if atom.get_parent().resname != 'HOH']
            self.natoms = len(self.atoms)
            #atom coordinates
            self.coordinates = np.array([atom.coord for atom in self.atoms])
            self.center = centroid(self.coordinates)
            self.coordinates -= self.center
            #atom color
            self.color = [colorrgba(restype(atom.get_parent().resname)) for atom in self.atoms]
            #atom radius
            self.radius = np.array([vrad(atom.get_id()) for atom in self.atoms])
            
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
            self.radius = np.array([vrad(atom.get_id()) for atom in self.atoms])
            
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
                self.n_atoms = len([atom for atom in self.residues[i] if atom.get_name() =='CA' or atom.get_name() == 'N'])
                self.color.append(np.tile(self.dsspcolor,(self.n_atoms,1)))
            if len(self.struct3)>1:
                self.color = np.concatenate(self.color)
            #atom radius
            self.radius = np.array([vrad(atom.get_id()) for atom in self.atoms])
            
VisPyViewer('data/1yd9.pdb',mode='cpk')