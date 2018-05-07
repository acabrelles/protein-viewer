from vispy import app, gloo, visuals, scene

from vispy.util.transforms import perspective, translate

import numpy as np

vertex = """

varying vec4  v_eye_position;
varying vec3  v_light_direction;

void main (void) {
    v_eye_position = $view * $model * vec4($position,1.0);
    v_light_direction = normalize($light_position);
    float dist = length(v_eye_position.xyz);

    gl_Position = $transform($projection * v_eye_position);

    vec4  proj_corner = $projection * vec4($radius, $radius, v_eye_position.z, v_eye_position.w);  // # noqa
    gl_PointSize = 512.0 * proj_corner.x / proj_corner.w;
}
"""

fragment = """

varying vec4  v_eye_position;
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
    pos.z += $radius*z;
    vec3 pos2 = pos.xyz;
    pos = $projection * pos;
    // gl_FragDepth = 0.5*(pos.z / pos.w)+0.5;
    vec3 normal = vec3(x,y,z);
    float diffuse = clamp(dot(normal, v_light_direction), 0.0, 1.0);

    // Specular lighting.
    vec3 M = pos2.xyz;
    vec3 O = v_eye_position.xyz;
    vec3 L = $light_spec_position;
    vec3 K = normalize(normalize(L - M) + normalize(O - M));
    // WARNING: abs() is necessary, otherwise weird bugs may appear with some
    // GPU drivers...
    float specular = clamp(pow(abs(dot(normal, K)), 40.), 0.0, 1.0);
    vec3 v_light = vec3(1., 1., 1.);
    gl_FragColor.rgb = (.15*$color + .55*diffuse * $color
                        + .35*specular * v_light);
}
"""

class MySpheresVisual(visuals.Visual):
    """Visual that draws a 3d plot
    
    Parameters
    ----------
    position: array of coordinates
    color: array of colors
    radius: array of radius
    W: width of canvas
    H: height of canvas
    """
    def __init__(self, coordinate, color, radius, W, H):
        visuals.Visual.__init__(self, vertex, fragment)
        
        self.size = W,H
        #Camera settings
        self.translate = 120
        self.view = translate((0,0, -self.translate), dtype=np.float32)
        self.model = np.eye(4, dtype=np.float32)
        self.projection = perspective(45.0, self.size[0] / float(self.size[1]), 1.0, 1000.0)
        
        self.shared_program.vert['model'] = self.model
        self.shared_program.vert['view'] = self.view
        self.shared_program.vert['projection'] = self.projection
        self.shared_program.frag['projection'] = self.projection
       
        self.shared_program.vert['light_position'] = 0., 0., 2.
        self.shared_program.frag['light_spec_position'] = -5., 5., -5.
        
        self.shared_program.vert['position'] = coordinate
        self.shared_program.vert['radius'] = radius
        self.shared_program.frag['radius'] = radius
        self.shared_program.frag['color'] = color
        
        self._draw_mode = 'points'
        
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

canvas = scene.SceneCanvas(keys='interactive', app='pyqt4', bgcolor='black',
                           size=(W, H), show=True)

view = canvas.central_widget.add_view()
view.camera = 'arcball'

spheres = [MySpheres(coordinates[0],color[0],radius[0],W,H, parent=view.scene)]

canvas.app.run()