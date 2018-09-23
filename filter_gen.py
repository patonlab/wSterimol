#!/usr/bin/python
from __future__ import print_function, absolute_import

import subprocess, sys, os, math, datetime, shutil
from os import listdir
from os.path import isfile, join

from pymol import cmd

########################################
########### F I L T E R - G ############
########################################

# Filter impossible structures after being generated via Pymol.
# Plausible structures can be visualized from visualize.py
# Use in Pymol command prompt:
# run log.py
# run sterimoltools.py
# run setup.py
# run filter_gen.py
# filter_gen (directory, setup_path, verbose)
# example: filter_gen conformers

def filter_gen(directory = "temp", setup_path = "default", verbose = "False"):
    # If the directory exists
    if os.path.exists(directory):
        # Log generation
        log_path = "log-%s" % (datetime.date.today())
        if os.path.exists(log_path+".pylog"):
            print("Warning: Continuing previous log file [%s.pylog]" % log_path)
        log = Log(log_path,"pylog")
        log.write("\n\n########################################\n########### F I L T E R - G ############\n########################################\n\n")
        #verbose
        if verbose.lower() in ['true', '1', 't', 'y', 'yes']: verbose = True
        else: verbose = False
        # Create a bin folder for the filtered structures
        bin_path = "%s/filter-bin-%s" % (directory, datetime.date.today())
        if not os.path.exists(bin_path):
            os.makedirs(bin_path)
            log.write("Filtration folder has been created in the directory to store the filtered structures [%s]" % bin_path, verbose)
        # Retrieve all the files in the directory
        files = [f for f in listdir(directory) if isfile(join(directory, f))]
        if len(files) > 0:
        # get setup.ini
            setup = Setup(log, setup_path)
            if setup.isLoaded() == True:
                #filter the files
                for file in files:
                    if accepted_file(file):
                        log.write("\nStart analysis file %s " % join(directory, file), verbose)
                        if nonbonded_clash(directory, file, setup.radii, setup.RJCT, log, verbose):
                            # There is a geometrical clash, the file is rejected
                            try:
                                shutil.move(join(directory, file), join(bin_path, file))
                            except:
                                log.write("Warning: Failed to put the file [%s] in the bin! Do it yourself!" % join(directory, file))
                                pass
                            log.write("*********** Rejected! Put in the bin!", verbose)
                        else: log.write("*********** Accepted!", verbose)
                log.write("----------------------------\n---- Normal Termination ----\n----------------------------\n")
                log.finalize()
            else: log.write("Error: Failed to load setup.ini in [%s]. Fix it to continue." % setup_path)
        else: log.write("Error: No file to filter in the directory [%s]" % directory)
    else: print("FATAL ERROR: Specified directory doesn't exist [%s]" % directory)


def nonbonded_clash(directory, file, radii, RJCT, log, verbose):
    # retrieve the data
    try: MolSpec = getinData(join(directory, file))
    except NameError: return False
    if hasattr(MolSpec, "ATOMTYPES"):
        if checkDists(MolSpec, radii, RJCT, log, verbose) != 0:
            return True
        else: return False
    else:
        log.write("Warning: Failed to read the file. No calculation. Skip the file. [%s]" % join(directory, file))
        return False

# Filter prior to optimization - if there are any very close nonbonded contacts, a non-zero value is returned
def checkDists(MolSpec, radii, RJCT, log, verbose):
 # RJCT is highly important. Define how close 2 atoms can be according to their diameter.
    checkval = 0
    for i in range(0,len(MolSpec.CARTESIANS)): # for each atom i
        bondedatomlist = [] # retrieve bonded atoms to i
        for partners in MolSpec.CONNECTIVITY[i]: bondedatomlist.append(int(partners.split("__")[0])-1)
        for j in range(i+1,len(MolSpec.CARTESIANS)): # check other atoms
            bond = 0
            for bondedatom in bondedatomlist:
                if j == bondedatom: bond = 1 # the j atom is bonded to atom i, ignore
            if bond == 0: # if not bonded, then check distance
                totdist = abs(calcdist(i, j, MolSpec.CARTESIANS))
                bump = RJCT*(atomicModelRadius(MolSpec.ATOMTYPES, MolSpec.CARTESIANS, radii, i)+atomicModelRadius(MolSpec.ATOMTYPES, MolSpec.CARTESIANS, radii, j))
                if totdist<bump: # distance is shorter than theoretical bumping
                    checkval = checkval+1
                    log.write("Rejecting structure! %3s %3s %3s %3s : distance = %6.3f Ang. bump = %6.3f Ang" % (MolSpec.ATOMTYPES[i],(i+1),MolSpec.ATOMTYPES[j],(j+1),totdist, bump), verbose)
    return checkval

# CPK VdW radii in pm
cpk_radii = [150,160,160,150,170,170,170,150,100,150,170,145,135,135,140,170,100,135,180,140,195,215]
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
                radius = cpk_radii[j]/100 # in angstroms
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

def calcdist(atoma,atomb,coords):
    x1=coords[atoma][0]
    y1=coords[atoma][1]
    z1=coords[atoma][2]
    x2=coords[atomb][0]
    y2=coords[atomb][1]
    z2=coords[atomb][2]
    ba = [x1-x2, y1-y2, z1-z2]
    dist = math.sqrt(ba[0]*ba[0]+ba[1]*ba[1]+ba[2]*ba[2])
    return dist

def atomicnumber(element):
    periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
                     "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
                     "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]
    
    for i in range(0,len(periodictable)):
        if element == periodictable[i]: return i

def accepted_file(file):
    filesplit = file.split(".")
    if len(filesplit) > 1: # there is an extension
        for suffix in ["pdb", "com"]: # authorized suffix
            if filesplit[len(filesplit)-1] == suffix: return True
    return False


cmd.extend("filter_gen",filter_gen)
