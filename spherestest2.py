from vispy import app, gloo, visuals, scene

from vispy.scene import ArcballCamera

from vispy.util.transforms import perspective, translate

import numpy as np

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

    gl_Position = $transform(u_projection * v_eye_position);
    //gl_Position = $transform(vec4(a_position, 1));

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

class MySpheresVisual(visuals.Visual):
    """Visual that draws a 3d plot
    
    Parameters
    ----------
    coordinates: array of coordinates
    color: array of colors
    radius: array of radius
    W: width of canvas
    H: height of canvas
    """
    def __init__(self, coordinates, color, radius, W, H):
        visuals.Visual.__init__(self, vertex, fragment)
        
        self.size = W,H
        #Camera settings
        self.translate = 120
        self.view = translate((0,0, -self.translate), dtype=np.float32)
        self.model = np.eye(4, dtype=np.float32)
        self.projection = np.eye(4, dtype=np.float32)
        
        self.apply_zoom()

        self.shared_program['u_model'] = self.model
        self.shared_program['u_view'] = self.view
        
        self.natoms = len(coordinates)
        
        #Loading data and type
        self.load_data()
        self._draw_mode = 'points'
    
    def apply_zoom(self):
        width, height = self.size        
        gloo.set_viewport(0, 0, width, height)
        self.projection = perspective(25.0, width / float(height), 50.0, 200.0)
        self.shared_program['u_projection'] = self.projection

    def load_data(self):
        
        """Make an array with all the data and load it into VisPy Gloo"""
        
        data = np.zeros(self.natoms, [('a_position', np.float32, 3),
                            ('a_color', np.float32, 4),
                            ('a_radius', np.float32, 1)])

        data['a_position'] = coordinates
        data['a_color'] = color
        data['a_radius'] = radius#*view.transforms.pixel_scale

        self.shared_program.bind(gloo.VertexBuffer(data))
        
        self.shared_program['u_light_position'] = 0., 0., 2.
        self.shared_program['u_light_spec_position'] = -5., 5., -5.
    
    def _prepare_transforms(self,view):
        view.view_program.vert['transform'] = view.get_transform()

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

MySpheres = scene.visuals.create_visual_node(MySpheresVisual)

canvas = scene.SceneCanvas(keys='interactive', app='pyqt4', bgcolor='white',
                           size=(W, H), show=True)

view = canvas.central_widget.add_view()

view.camera = ArcballCamera()#(fov=20, distance=300)

spheres = [MySpheres(coordinates,color,radius,W,H,parent=view.scene)]

canvas.app.run()