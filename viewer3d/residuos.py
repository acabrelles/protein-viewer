from direct.showbase.ShowBase import ShowBase
from panda3d.core import ColorAttrib
from panda3d.core import LVecBase4f
from panda3d.core import NodePath, PandaNode
from panda3d.core import AmbientLight, DirectionalLight

from Bio.PDB import *

from molecular_data import resdict, colorrgba, restype, vrad


#Creamos la base de la ventana
base = ShowBase()

#Desglosamos archivo PDB
parser = PDBParser(QUIET=True,PERMISSIVE=True)
structure = parser.get_structure('MACROH2A', 'data/1yd9.pdb')

#Creamos el modelo

pnode = render.attachNewNode("Model")

#Listado de residuos, solo contando los aminoacidos
residues = [residue for residue in structure.get_residues() if residue.get_resname() in resdict.keys()]
for residue in residues:
    resid = residue.get_resname()
    color=colorrgba(restype(resid))
    atoms = [atom for atom in residue.get_atoms()]
    for atom in atoms:
        x,y,z=atom.coord
        id=atom.get_id()
        a = loader.loadModel("data/atom_sphere")
        a.setPos(x, y, z)
        a.setColor(color)
        a.setScale(vrad(id))
        a.reparentTo(pnode)

pnode.flattenStrong

#Seleccionamos el resto de 'residuos' que ha determinado el parser, sin incluir aguas
residues2 = [residue for residue in structure.get_residues() if not residue in residues and residue.get_resname() != 'HOH']
for residue in residues2:
    atoms = [atom for atom in residue.get_atoms()]
    for atom in atoms:
        x,y,z=atom.coord
        id=atom.get_id()
        a = loader.loadModel("data/atom_sphere")
        a.setPos(x, y, z)
        a.setColor(colorrgba(id))
        a.setScale(vrad(id))
        a.reparentTo(pnode)

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
