#!/usr/bin/python

import subprocess, sys, os, math, datetime, shutil, getpass

from pymol import cmd

########################################
########## V I S U A L I Z E ###########
########################################

# Create functions to visualize the different created conformers in Pymol
# Use in Pymol command prompt:
# run visualize.py

# AddConformers
# example: AddConformers
def AddConformers(path = "temp"):
    if os.path.exists(path):
        dir_files = [file for file in os.listdir(path) if os.path.isfile(os.path.join(path, file))]  # existing files
        if len(dir_files) > 0:
            for file_name in dir_files:
                filesplit = file_name.split(".")
                if len(filesplit) > 1: # there is an extension
                    if filesplit[-1] == "pdb": # only pdb is accepted here
                        cmd.load(os.path.join(path, file_name)) #load the structures
                        BallnStick( filesplit[0] ) # make it look pretty
        else: print "Error: No files to load in directory [%s]" % path
    else: print "Error: The path doesn't exist [%s]" % path

cmd.extend("AddConformers", AddConformers)

# CPK van der Waals radii in Angstrom
cpk_radii = [1.50,1.60,1.60,1.50,1.70,1.70,1.70,1.50,1.00,1.50,1.70,1.45,1.35,1.35,1.40,1.70,1.00,1.35,1.80,1.40,1.95,2.15]
# Verloop's original Sterimol parameters use CPK atomic VdW radii based on atom-type definitions
sterimol_atomtypes = ["C", "C2", "C3", "C4", "C5/N5", "C6/N6", "C7", "C8", "H", "N", "C66", "N4", "O", "O2", "P", "S", "S1", "F", "C1", "S4", "B1", "I"]
## covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197 ##
## values for metals decreased by 10 % ##
rcov = [0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,
        1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54,
        1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09,
        1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39,
        1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26,
        1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57,
        1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53,
        1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32,
        1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58,
        1.52, 1.53, 1.54, 1.55]
#Bondi van der Waals radii for all atoms from: Bondi, A. J. Phys. Chem. 1964, 68, 441-452, except hydrogen, which is taken from Rowland, R. S.; Taylor, R. J. Phys. Chem. 1996, 100, 7384-7391
#Radii that are not available in either of these publications have RvdW = 2.00 Angstrom
bondi_radii = [0.0,1.09, 1.40, 1.82,2.00,2.00,1.70,1.55,1.52,1.47,1.54,2.27,1.73,2.00,2.10,1.80,1.80,1.75,1.88,2.75,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.40,1.39,1.87,2.00,1.85,1.90,
        1.85,2.02,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.72,1.58,1.93,2.17,2.00,2.06,1.98,2.16,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.72,1.66,1.55,1.96,2.02,2.00,2.00,2.00,
        2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.86]
#Some useful arrays for chemists
periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr", "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl", "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo","Bq"]

    
def atomicModelRadius(ATOMTYPES, CARTESIANS , radii, atomid):
    if radii == "cpk": # CPK needs to know the "atomtype", and the atomtype is defined according to the atom environment... Not very efficient
        atomic_co_no = ncoord_atom(atomid, ATOMTYPES, CARTESIANS)
        sterimol_type = generate_atom_types_atom(ATOMTYPES[atomid], atomic_co_no)
        # find the right atomtypes corresponding to our atom
        radius = 2.0 # doesn't recognise it, 2.0 angstroms by default
        for j in range(0,len(sterimol_atomtypes)):
            if sterimol_type == sterimol_atomtypes[j]:
                radius = cpk_radii[j]
                break
    elif radii == "bondi":
        massno = atomicnumber(ATOMTYPES[atomid])
        if massno<len(bondi_radii): radius = bondi_radii[massno]
        else: radius = 2.0
    return radius

# Generate Sterimol atom type from connectivity data
def generate_atom_types_atom(atom, cn):
    if atom == "H": return "H"
    elif atom == "P": return "P"
    elif atom == "F": return "F"
    elif atom == "Cl": return "C1"
    elif atom == "Br": return "B1"
    elif atom == "I": return "I"
    elif atom == "O": #Sterimol distinguishes between "normal", and double-bonded O atoms
        if cn < 1.5: return "O2"
        if cn > 1.5: return "O"
    elif atom == "S": #Sterimol distinguishes between "normal", tetrahedral, and octohedral S atoms
        if cn < 2.5: return "S"
        if 5.5 > cn > 2.5: return "S4"
        if cn > 5.5: return "S1"
    elif atom == "N": #Sterimol distinguishes between tetrahedral and planar (amide) N atoms
        if cn > 2.5: return "N"
        if cn < 2.5: return "C6/N6"
    elif atom == "C": #Sterimol distinguishes between myriad types of C atoms ...
        if cn < 2.5: return "C3"
        if 3.5 > cn > 2.5: # need to differentiate between sp2 carbon and aromatic carbon ...
            return "C6/N6" # assumes aromatic rather than sp2
        if cn > 3.5: return "C"

# Calculation of atomic coordination numbers (taken from Grimme's DFTD3 definitions) for 1 atom
def ncoord_atom(atomid, atomtype, coords):
    max_elem = 94
    k1 = 16.0
    k2 = 4.0/3.0
    xn = 0.0
    for iat in range(0,len(coords)):
        if iat != atomid:
            dx = coords[iat][0] - coords[atomid][0]
            dy = coords[iat][1] - coords[atomid][1]
            dz = coords[iat][2] - coords[atomid][2]
            r2 = dx*dx+dy*dy+dz*dz
            r = math.pow(r2,0.5)
            Zi = -1
            Ziat = -1
            for k in range(0,max_elem):
                if atomtype[atomid].find(elements[k])>-1:Zi=k
                if atomtype[iat].find(elements[k])>-1:Ziat=k
                if Zi > -1 and Ziat > -1: break
            rco = rcov[Zi]+rcov[Ziat]
            rco = rco*k2
            rr=rco/r
            damp=1.0/(1.0+math.exp(-k1*(rr-1.0)))
            xn=xn+damp
    return xn

elements = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si",
        "P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni",
        "Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo",
        "Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba",
        "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
        "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
        "Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
        "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq",
        "Uup","Uuh","Uus","Uuo"]
    
def atomicnumber(element):
    for i in range(1,len(periodictable)):
        if element == periodictable[i]: return i
    
   
def BallnStick( arg1 ):
    cmd.show("sticks", arg1)
    cmd.show("spheres", arg1)
    cmd.color("gray","elem C and "+arg1)
    cmd.set("stick_radius",0.07, arg1)
    cmd.set("sphere_scale",0.15, arg1)
    cmd.alter(arg1+" and elem H", "vdw=0.75")
    cmd.set("stick_color","black", arg1)
    cmd.set("dash_gap",0.01, arg1)
    cmd.hide("nonbonded", arg1)
    cmd.hide("lines", arg1)

cmd.extend( "BallnStick", BallnStick );

def Add_VDW( object, radii = ""):
    cmd.copy(object+"_vdw", object)
    cmd.alter(object+"_vdw and elem H", "vdw=1.09")
    radii = radii.lower()
    if radii == "cpk" or radii == "bondi":
        # retrieve the data (coords and atomtypes)
        ATOMTYPES = []
        CARTESIANS = []
        atomslist = cmd.get_model(object+"_vdw")
        for i in range(0,len(atomslist.atom)):
            ATOMTYPES.append(atomslist.atom[i].symbol)
            CARTESIANS.append(atomslist.atom[i].coord)
        for i in range(0,len(ATOMTYPES)):
            atomid = i + 1 # careful, pymol numbering starts at 1, but list starts at 0
            cmd.alter("%s_vdw & id %s" % (object, atomid), "vdw=%.2f" % atomicModelRadius(ATOMTYPES, CARTESIANS , radii, i))
    cmd.rebuild()
    cmd.set("sphere_scale", 1, object+"_vdw")
    cmd.hide("nonbonded", object+"_vdw")
    cmd.hide("lines", object+"_vdw")
    cmd.hide("sticks", object+"_vdw")
    cmd.set("sphere_transparency", 0.7, object+"_vdw")
 
cmd.extend( "Add_VDW", Add_VDW );

# Setup a nice environment for the Paton group images
cmd.bg_color("white")
cmd.set("ray_opaque_background", "off")
cmd.set("specular", 0.25)
cmd.set("spec_power", 300)
cmd.set("spec_reflect", 0.5)
cmd.util.ray_shadows("light")

cmd.set("antialias", 1)
cmd.set("orthoscopic", 0)
cmd.set("field_of_view", 30)
cmd.set("transparency", 0.5)






