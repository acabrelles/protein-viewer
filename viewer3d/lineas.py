from direct.showbase.ShowBase import ShowBase
from panda3d.core import ColorAttrib
from panda3d.core import LVecBase4f
from panda3d.core import NodePath, PandaNode
from panda3d.core import AmbientLight, DirectionalLight
from panda3d.core import LineSegs

from Bio.PDB import *

import numpy as np

from molecular_data import colorrgba, vrad


#Creamos la base de la ventana
base = ShowBase()

#Desglosamos archivo PDB
parser = PDBParser(QUIET=True,PERMISSIVE=True)
structure = parser.get_structure('MACROH2A', 'data/1yd9.pdb')

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
		a = loader.loadModel("data/atom_sphere")
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
