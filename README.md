# protein-viewer
Using matplotlib, VisPy and Panda3D for viewing proteins through their PDB archives.

This project has as an objective using different libraries for viewing 3D protein models.

The libraries are:

* Mathematical library, matplotlib
* OpenGL based library, VisPy
* Videogames engine library, Panda3D. This one is extra to the project.

### Requirements
* Numpy
* BioPython
* VisPy and PyQt4
* Tkinter
* Panda3D (extra)

# Commands
There's 4 visualizations per library:

* CPK: Atoms are coloured by the CPK standard.
* Aminoacid: Atoms change color depending on the type of residue they belong to.
* Backbone: Chains representations, showing only the CA-N bonds of the skeleton.
* DSSP: Protein structure prediction, showing only the CA-N bonds of the skeleton.

In the data folder there's a PDB file provided for testing, `1yd9.pdb`

## Tkinter
To load the GUI, execute the next command:

```
$ python pviewer.py
```

A window with the settings to fiddle with will appear.

## MatViewer
When importing this package, you can use the next command to represent your protein:

```
MatViewer(data.pdb, [cpk|aminoacid|backbone|dssp])
``` 

## VisPyViewer
When importing this package, you can use the next command to repressent your protein:

```
VisPyViewer(data.pdb, [cpk|aminoacid|backbone|dssp])
```

## PANDA3D - WIP
Execute the command for the example: 

```
$ python pandaviewer.py data/1yd9.pdb
```

Controls are:

1: CPK

2: Aminoacid type

3: Backbone

4: DSSP

c: chain cloud

Left arrow key: Horizontal rotation

Up arrow key: Vertical rotation

Down arrow key: Stop rotation


(Camera is a bit wonky)

Left click for panning

Middle click for rotation 1

Right click for zooming

Middle + Right click for rotation 2

Escape: Exit
