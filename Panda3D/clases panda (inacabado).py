from direct.showbase.ShowBase import ShowBase
from panda3d.core import ColorAttrib
from panda3d.core import LVecBase4f
from panda3d.core import NodePath, PandaNode
from panda3d.core import AmbientLight, DirectionalLight
from panda3d.core import LineSegs

from Bio.PDB import *

import dict

import numpy as np

class panda3d(ShowBase):

	modes={'CPK':'load_CPK', 'backbone':'load_bbone'}

	def __init(self, pdb_file, mode='CPK'):
		ShowBase.__init__(self)
		if mode not in modes.keys():
			raise Exception('No se reconoce el modo de visualizacion %s' % mode)
		self.mode = mode

		#Desglosamos archivo PDB
		parser = PDBParser(QUIET=True,PERMISSIVE=True)
		structure = parser.get_structure('Structure', pdb_file)	
	
		#Creamos el modelo
		self.pnode = render.attachNewNode("Model")
		if mode == 'CPK':
			self.load_CPK(structure, self.pnode)
		elif mode == 'backbone':
			self.load_bbone(structure, self.pnode)
		#elif mode == 'elquesea':

		#Colocamos la proteina en el centro

		self.pnode.setPos(0,0,0)

		#Colocamos la camara en el centro
				
		p_radius= self.pnode.getBounds().getRadius()
		p_center= self.pnode.getBounds().getCenter()		
		xc,yc,zc = p_center
		
		base.cam.setPos(xc,-150-yc-2*p_radius,zc)
		base.cam.lookAt(xc,yc,zc)

		#Luz ambiental
		
		self.ambiente = AmbientLight('aluz')
		self.ambiente.setColor(LVecBase4f(0.16, 0.16, 0.17, 1.0))
		self.luza = render.attachNewNode(ambiente)
		render.setLight(self.luza)

		#Luz direccional
		
		self.direccional = DirectionalLight('dluz')
		self.direccional.setColor(LVecBase4f(0.8,0.7,0.75,1.0))
		self.direccional.setShadowCaster(True,512,512)
		render.setShaderAuto()
		self.luzd = render.attachNewNode(direccional)
		self.luzd.setPos(0,-50,0)
		render.setLight(self.luzd)
		self.luzd.lookAt(xc,yc,zc)

	#Definimos modo de representacion
	def load_CPK(self, structure, node):
		for atom in structure.get_atoms():
			x,y,z=atom.coord
			id=atom.get_id()
			a = loader.loadModel("atom_sphere")
			a.setPos(x, y, z)
			a.setColor(dict.colorrgba(id))
			a.setScale(dict.vrad(id))
			a.reparentTo(node)

		node.flattenStrong()

	def load_bbone(self, structure, node):
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
				a.setScale(dict.vrad(id)/2.5)

			lines = LineSegs()
			lines.setColor(ccolor)
			lines.moveTo(can_coordinates[0][0],can_coordinates[0][1],can_coordinates[0][2])
			for i in range(len(can_atoms))[1:]:
				lines.drawTo(can_coordinates[i][0],can_coordinates[i][1],can_coordinates[i][2])
			lines.setThickness(6)
			lnode = lines.create()
			linenp = NodePath(lnode)
			linenp.reparentTo(pnode)

		node.flattenStrong()