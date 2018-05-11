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


Spheres = scene.visuals.create_visual_node(SpheresVisual)


if __name__ == '__main__':
    from vispy import scene

    canvas = scene.SceneCanvas(keys='interactive', app='pyqt4', bgcolor='white',
                            size=(1000, 800), show=True)
    view = canvas.central_widget.add_view()
    view.camera = scene.ArcballCamera(fov=70, distance=50)

    N = 2000
    coordinates = np.random.normal(size=(N, 3), scale=10)
    color = np.random.random(size=(N, 4))
    color[:,3] = 1.0
    radius = np.random.normal(size=N, loc=4.0, scale=0.5)

    spheres = [Spheres(coordinates,color,radius,parent=view.scene)]

    #canvas.app.run()