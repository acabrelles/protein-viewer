from direct.showbase.ShowBase import ShowBase
from panda3d.core import ColorAttrib
from panda3d.core import LVecBase4f
from panda3d.core import NodePath, PandaNode
from panda3d.core import AmbientLight, DirectionalLight
from panda3d.core import LineSegs

from Bio.PDB import *

import numpy as np


#Creamos la base de la ventana
base = ShowBase()

#Desglosamos archivo PDB
parser = PDBParser(QUIET=True,PERMISSIVE=True)
structure = parser.get_structure('MACROH2A', '1yd9.pdb')

#Funcion para sacar radios

def vrad(atom):
    radius = {'H':1.2,
          'C':1.7, 'CA':1.7,'CB':1.7,'CG':1.7,'CZ':1.7,
          'N':1.55,
          'O':1.52,
          'F':1.47,'Cl':1.75,
          'Br':1.85,
          'I':1.98,
          'He':1.40, 'Ne':1.54, 'Ar':1.88,'Xe':2.16,'Kr':2.02,
          'P':1.8,
          'S':1.8,
          'B':1.92,
          'Li':1.82, 'Na':2.27,'K':2.75,'Rb':3.03,'Cs':3.43,
          'Be':1.53,'Mg':1.73,'Ca':2.31,'Sr':2.49,'Ba':2.68,'Ra':2.83}
    try:
        return radius[atom]
    except:
        return 1.5

#Funciones para sacar colores

def cpkcolor(atom):
    colors = {'H':'white',
          'C':'black', 'CA':'black','CB':'black','CG':'black','CZ':'black',
          'N':'blue',
          'O':'red',
          'F':'green','Cl':'green',
          'Br':'brown',
          'I':'darkviolet',
          'He':'turquoise', 'Ne':'turquoise', 'Ar':'turquoise','Xe':'turquoise','Kr':'turquoise',
          'P':'orange',
          'S':'yellow',
          'B':'salmon',
          'Li':'purple','Na':'purple','K':'purple','Rb':'purple','Cs':'purple',
          'Be':'darkgreen','Mg':'darkgreen','Ca':'darkgreen','Sr':'darkgreen','Ba':'darkgreen','Ra':'darkgreen',
          'Ti':'grey',
          'Fe':'darkorange'}
    try:
        return colors[atom]
    except:
        return 'pink'

def crgba(color):
    rgba = {'white':(1.0,1.0,1.0,1.0), 
        'black':(0.05,0.05,0.05,1.0),
        'blue':(0.0,0.0,1.0,1.0),
        'red':(1.0,0.0,0.0,1.0),
        'green':(0.13,0.78,0.0,1.0),
        'brown':(0.4,0.0,0.0,1.0),
        'darkviolet':(0.24,0.0,0.4,1.0),
        'turquoise':(0.0,0.78,0.84,1.0),
        'orange':(0.84,0.53,0.0,1.0),
        'yellow':(0.86,0.9,0.0,1.0),
        'salmon':(1.0,0.75,0.51,1.0),
        'purple':(0.35,0.0,0.59,1.0),
        'darkgreen':(0.0,0.35,0.0,1.0),
        'grey':(0.59,0.59,0.59,1.0),
        'darkorange':(0.86,0.45,0.0,1.0),
        'pink':(0.94,0.55,1.0,1.0)}
    try:
        return rgba[color]
    except:
        return rgba['pink']

#Creamos el modelo

pnode = render.attachNewNode("Model")

for chain in structure.get_chains():
	carr = np.random.rand(3,1)
	ccolor = float(carr[0]),float(carr[1]),float(carr[2]),1.0
	can_atoms = [atom for atom in chain.get_atoms() if atom.get_name() == 'CA' or atom.get_name() == 'N']
	can_coordinates = [atom.coord for atom in can_atoms]
	for atom in can_atoms:
		x,y,z = atom.coord
		id=atom.get_id()
		a = loader.loadModel("atom_sphere")
		a.setPos(x,y,z)
		a.reparentTo(pnode)
		a.setColor(ccolor)
		a.setScale(vrad(id)/2.5)

	lines = LineSegs()
	lines.setColor(ccolor)
	lines.moveTo(can_coordinates[0][0],can_coordinates[0][1],can_coordinates[0][2])
	for i in range(len(can_atoms))[1:]:
		lines.drawTo(can_coordinates[i][0],can_coordinates[i][1],can_coordinates[i][2])
	lines.setThickness(6)
	lnode = lines.create()
	linenp = NodePath(lnode)
	linenp.reparentTo(pnode)

pnode.flattenStrong()

#Colocamos la proteina en el centro

pnode.setPos(0,0,0)

#Colocamos la camara en el centro

p_radius= pnode.getBounds().getRadius()
p_center= pnode.getBounds().getCenter()
xc,yc,zc = p_center

base.cam.setPos(xc,-150-yc-2*p_radius,zc)
base.cam.lookAt(xc,yc,zc)

#Creamos iluminacion del ambiente (para la sombra)

ambiente = AmbientLight('aluz')
ambiente.setColor(LVecBase4f(0.16, 0.16, 0.17, 1.0))
luza = render.attachNewNode(ambiente)
render.setLight(luza)

#Creamos una iluminacion direccional (para dar volumen)

direccional = DirectionalLight('dluz')
direccional.setColor(LVecBase4f(0.8,0.7,0.75,1.0))
direccional.setShadowCaster(True,512,512)
render.setShaderAuto()
luzd = render.attachNewNode(direccional)
luzd.setPos(0,-50,0)
render.setLight(luzd)
luzd.lookAt(xc,yc,zc)

base.run()