resdict={'SER':'np','THR':'np','GLN':'np','ASN':'np','TYR':'np','CYS':'np',
         'ALA':'nnp','VAL':'nnp','LEU':'nnp','ILE':'nnp','MET':'nnp','PRO':'nnp','PHE':'nnp','TRP':'nnp','GLY':'nnp',
         'ASP':'neg','GLU':'neg',
         'LYS':'pos','ARG':'pos','HIS':'pos'}

colors={'H':'white',
        'C':'black', 'CA':'black','CB':'black','CG':'black','CZ':'black',
        'N':'blue',
        'O':'red',
        'F':'green','CL':'green',
        'BR':'brown',
        'I':'darkviolet',
        'HE':'turquoise', 'NE':'turquoise', 'AR':'turquoise','XE':'turquoise','KR':'turquoise',
        'P':'orange',
        'S':'yellow',
        'B':'salmon',
        'LI':'purple','NA':'purple','K':'purple','RB':'purple','CS':'purple',
        'BE':'darkgreen','MG':'darkgreen','Ca':'darkgreen','SR':'darkgreen','BA':'darkgreen','RA':'darkgreen',
        'TI':'grey',
        'FE':'darkorange',
        'np':'rebeccapurple',
        'nnp':'gray',
        'neg':'navy',
        'pos':'darkred'}

rgba = {'white':(1.0,1.0,1.0,1.0), 
	'black':(0.02,0.02,0.02,1.0),
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
	'pink':(0.94,0.55,1.0,1.0),
	'rebeccapurple':(0.34,0.25,0.63,1.0),
        'gray':(0.75,0.75,0.75,1.0),
	'navy':(0.0,0.06,0.51,1.0),
	'darkred':(0.55,0.0,0.0,1.0)}

#Sacados de https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
radius = {'H':1.2,
          'C':1.7, 'CA':1.7,'CB':1.7,'CG':1.7,'CZ':1.7,
          'N':1.55,
          'O':1.52,
          'F':1.47,'CL':1.75,
          'BR':1.85,
          'I':1.98,
          'HE':1.40, 'NE':1.54, 'AR':1.88,'XE':2.16,'KR':2.02,
          'P':1.8,
          'S':1.8,
          'B':1.92,
          'LI':1.82, 'NA':2.27,'K':2.75,'RB':3.03,'CS':3.43,
          'BE':1.53,'MG':1.73,'CA':2.31,'SR':2.49,'BA':2.68,'RA':2.83}

def restype(residue):
    return resdict[residue]

def color(atom):
    try:
        return colors[atom]
    except:
        return 'pink'

def crgba(color):
    return rgba[color]

def colorrgba(atom):
    try:
        return rgba[colors[atom]]
    except:
        return rgba['pink']

def vrad(atom):
    try:
        return radius[atom]
    except:
        return 1.5