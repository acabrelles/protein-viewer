from direct.showbase.ShowBase import ShowBase
from panda3d.core import ColorAttrib, TransparencyAttrib, AntialiasAttrib
from panda3d.core import LVecBase4f
from panda3d.core import NodePath, PandaNode, TextNode
from panda3d.core import AmbientLight, DirectionalLight
from direct.task import Task
from panda3d.core import LineSegs
from panda3d.core import Mat4

from Bio.PDB import PDBParser,DSSP

from math import pi, sin, cos
import numpy as np

from molecular_data import resdict, restype, colorrgba, vrad, crgbaDSSP

import sys, os

class PandaViewer(ShowBase):
    
    if len(sys.argv[1:])!=1:
        raise SystemExit("Uso: %s Archivo_PDB" % os.path.basename(sys.argv[0]))
    
    def __init__(self):
        ShowBase.__init__(self)
        self.cloud = False
        self.help = False
        self.screen_text = []
        
        #Desglosamos archivo PDB
        pdbdata = sys.argv[1]
        parser = PDBParser(QUIET=True,PERMISSIVE=True)
        structure = parser.get_structure('model', pdbdata)
        
        #Hacemos la prediccion DSSP
        model = structure[0]
        dssp = DSSP(model, pdbdata)
        
        #Creamos los modelos
        self.cpknode = render.attachNewNode("CPK")
        self.aanode = render.attachNewNode("Aminoacids")
        self.bbnode = render.attachNewNode("BackBone")
        self.dsspnode = render.attachNewNode("DSSP")
        self.nnode = render.attachNewNode("Cloud")
        
        #CPK
        for atom in structure.get_atoms():
            x, y, z = atom.coord
            atomid = atom.get_id()
            a = loader.loadModel("data/atom_sphere")
            a.setPos(x, y, z)
            a.reparentTo(self.cpknode)
            a.setColor(colorrgba(atomid))
            a.setScale(vrad(atomid))

        self.cpknode.flattenStrong()
            
        #Aminoacids
        self.residues = [residue for residue in structure.get_residues() if residue.get_resname() in resdict.keys()]
        for residue in self.residues:
            resid = residue.get_resname()
            color = colorrgba(restype(resid))
            atoms = [atom for atom in residue.get_atoms()]
            for atom in atoms:
                x, y, z=atom.coord
                atomid=atom.get_id()
                a = loader.loadModel("data/atom_sphere")
                a.setPos(x, y, z)
                a.setColor(color)
                a.setScale(vrad(atomid))
                a.reparentTo(self.aanode)

        self.residues2 = [residue for residue in structure.get_residues() if not residue in self.residues and residue.get_resname() != 'HOH']
        for residue in self.residues2:
            atoms = [atom for atom in residue.get_atoms()]
            for atom in atoms:
                x, y, z=atom.coord
                atomid=atom.get_id()
                a = loader.loadModel("data/atom_sphere")
                a.setPos(x, y, z)
                a.setColor(colorrgba(atomid))
                a.setScale(vrad(atomid))
                a.reparentTo(self.aanode)
        self.aanode.flattenStrong()
        self.aanode.hide()
        
        #Backbone
        for chain in structure.get_chains():
            carr = np.random.rand(3,1)
            ccolor = float(carr[0]),float(carr[1]),float(carr[2]),1.0
            can_atoms = [atom for atom in chain.get_atoms() if atom.get_name() == 'CA' or atom.get_name() == 'N']
            can_coordinates = [atom.coord for atom in can_atoms]
            for atom in can_atoms:
                x, y, z = atom.coord
                atomid=atom.get_id()
                a = loader.loadModel("data/atom_sphere")
                a.setPos(x,y,z)
                a.reparentTo(self.bbnode)
                a.setColor(ccolor)
                a.setScale(vrad(atomid)/2.5)

            lines = LineSegs()
            lines.setColor(ccolor)
            lines.moveTo(can_coordinates[0][0],can_coordinates[0][1],can_coordinates[0][2])
            for i in range(len(can_atoms))[1:]:
                lines.drawTo(can_coordinates[i][0],can_coordinates[i][1],can_coordinates[i][2])
            lines.setThickness(6)
            lnode = lines.create()
            self.linenp = NodePath(lnode)
            self.linenp.instanceTo(self.bbnode)

            #Cloud
            catoms = [atom for atom in chain.get_atoms()]
            for atom in catoms:
                x, y, z = atom.coord
                atomid=atom.get_id()
                a = loader.loadModel("data/atom_sphere")
                a.setPos(x,y,z)
                a.reparentTo(self.nnode)
                a.setColor(ccolor)
                a.setScale(vrad(atomid)*1.1)

        self.bbnode.flattenStrong()
        self.bbnode.hide()
        self.nnode.setTransparency(TransparencyAttrib.MAlpha)
        self.nnode.setAlphaScale(0.3)
        self.nnode.hide()
        
        #DSSP
        self.linenp.instanceTo(self.dsspnode)
        self.struct3 = [dssp[key][2] for key in list(dssp.keys())]    
    
        for i in range(len(self.struct3)):
            dsspcolor = crgbaDSSP(self.struct3[i])
            can_atoms = [atom for atom in self.residues[i] if atom.get_name() == 'CA' or atom.get_name() == 'N']
            for atom in can_atoms:
                x, y, z = atom.coord
                atomid=atom.get_id()
                a = loader.loadModel("data/atom_sphere")
                a.setPos(x, y, z)
                a.reparentTo(self.dsspnode)
                a.setColor(dsspcolor)
                a.setScale(vrad(atomid)/2.5)
            self.dsspnode.flattenStrong()
            self.dsspnode.hide()
        
        #Colocamos la proteina en el centro
        self.cpknode.setPos(0,0,0)
        self.bbnode.setPos(0,0,0)
        self.aanode.setPos(0,0,0)
        self.nnode.setPos(0,0,0)
        
        #Colocamos la camara en el centro
        xc, yc, zc = self.cpknode.getBounds().getCenter()
        self.center = xc, yc, zc
        self.pradius = self.cpknode.getBounds().getRadius()
        self.center_camera()
        
        #Creamos la iluminacion de ambiente
        self.ambient = AmbientLight('alight')
        self.ambient.setColor(LVecBase4f(0.16, 0.16, 0.17, 1.0))
        self.alight = render.attachNewNode(self.ambient)
        render.setLight(self.alight)
        
        #Creamos la iluminacion direccional
        self.directional = DirectionalLight('dlight')
        self.directional.setColor(LVecBase4f(0.8, 0.7, 0.75, 1.0))
        self.directional.setShadowCaster(True,512,512)
        render.setShaderAuto()
        self.dlight = render.attachNewNode(self.directional)
        self.dlight.setPos(0,-50,0)
        render.setLight(self.dlight)
        self.dlight.lookAt(self.cpknode.getBounds().getCenter())
        
        # Post procesado      
        render.setAntialias(AntialiasAttrib.MAuto)
        
        #Teclado
        self.accept('c', self.toggle_cloud)
        self.accept('1', self.showmodel, [self.cpknode])
        self.accept('2', self.showmodel, [self.aanode])
        self.accept('3', self.showmodel, [self.bbnode])
        self.accept('4', self.showmodel, [self.dsspnode])
        self.accept('x', self.center_camera)
        self.accept('arrow_left', self.taskMgr.add, [self.spinCameraTaskX, "SpinCameraTaskX"])
        self.accept('arrow_up', self.taskMgr.add, [self.spinCameraTaskY, "SpinCameraTaskY"])
        self.accept('arrow_down', self.stop_camera)
        self.accept('escape', sys.exit)
        
    def center_camera(self):
        base.cam.setPos(self.center[0], -10-self.center[1]-4*self.pradius, self.center[2])
        base.cam.lookAt(self.center)
        
    def toggle_cloud(self):
        self.cloud = not self.cloud
        if self.cloud:
            self.nnode.show()
        else:
            self.nnode.hide()
            
    def spinCameraTaskX(self,task):
        base.disableMouse()
        taskMgr.remove("SpinCameraTaskY")
        angleDegrees = task.time * 20.0
        angleRadians = angleDegrees * (pi / 180.0)
        self.camera.lookAt(self.center)
        self.camera.setHpr(angleDegrees, 0, 0)
        return Task.cont

    def spinCameraTaskY(self,task):
        base.disableMouse()
        taskMgr.remove("SpinCameraTaskX")
        angleDegrees = task.time * 20.0
        angleRadians = angleDegrees * (pi / 180.0)
        self.camera.lookAt(self.center)
        self.camera.setHpr(0, angleDegrees, 0)
        return Task.cont
    
    def stop_camera(self):
        mat=Mat4(camera.getMat())
        mat.invertInPlace()
        base.mouseInterfaceNode.setMat(mat)
        base.enableMouse()
        self.taskMgr.remove("SpinCameraTaskX")
        self.taskMgr.remove("SpinCameraTaskY")
        
    def showmodel(self, node):
        self.cpknode.hide()
        self.aanode.hide()
        self.bbnode.hide()
        self.dsspnode.hide()
        node.show()
            
        
app = PandaViewer()
app.run()