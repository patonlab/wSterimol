#!/usr/bin/python


###Sterimol (and Tolman CA) Calculator###

###############################################################
#                       sterimoltools.py                      #
#                                                             #
###############################################################


#Python Libraries
import subprocess, sys, os
from numpy import *
from scipy import *
from math import *
import numpy as np
#from vpython import *

#Chemistry Libaries
#from radialdata import *
#from pars import *

#Avoid number error warnings
import warnings
warnings.filterwarnings("ignore")

#Chemistry Arrays
periodictable = ["Bq","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
             "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
             "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

metals = ["Li","Be","Na","Mg","Al","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Uut","Fl","Uup","Lv"]

# Verloop's original Sterimol parameters use CPK atomic VdW radii based on atom-type definitions
sterimol_atomtypes = ["C", "C2", "C3", "C4", "C5/N5", "C6/N6", "C7", "C8", "H", "N", "C66", "N4", "O", "O2", "P", "S", "S1", "F", "C1", "S4", "B1", "I"]

# CPK VdW radii in pm
cpk_radii = [150,160,160,150,170,170,170,150,100,150,170,145,135,135,140,170,100,135,180,140,195,215]

def getfragment(atom,molcart):
   bondlist=[atom]
   for a in range(len(molcart)):
      if calcdist(atom,a,molcart)<1.92 and a not in bondlist:bondlist.append(a)

      for b in range(len(bondlist)):
         for c in range(len(molcart)):

            if calcdist(bondlist[b],c,molcart)<1.92 and c not in bondlist:bondlist.append(c)
   return bondlist

def connectivity(atom,molcart,aty):
   con=[]
   for a in range(len(molcart)):
      if aty[a]in metals and molcart[a] != molcart[atom] and 0.1<calcdist(a,atom,molcart)<2:con.append(a)
      if molcart[a] != molcart[atom] and 0.1<calcdist(a,atom,molcart)<1.7:con.append(a)
   return len(con)
def genradii(atom,molcart,aty):
   #molcart=fileData.CARTESIANS
   con=connectivity(atom,molcart,aty)
   if con==0:con=1
   type=aty[atom]
   arow=periodictable.index(type)
   radius=molmod[arow][con]
   if radius==0:radius=1;print "Warning: No atomic radii found", arow, con
   return radius

def rotrel(vect1,vect2,vect3):
   ax=np.cross(vect1,vect2)
   ang=math.acos((np.dot(vect1,vect2))/(np.linalg.norm(vect1)*np.linalg.norm(vect2)))
   norm=1/(np.linalg.norm(ax))
   axnorm=np.dot(ax,norm)
   ux=axnorm[0]
   uy=axnorm[1]
   uz=axnorm[2]
   a=math.cos(ang)+((ux*ux)*(1-math.cos(ang)))
   b=(ux*uy*(1-math.cos(ang)))-(uz*math.sin(ang))
   c=(ux*uz*(1-math.cos(ang)))+(uy*math.sin(ang))
   d=(uy*ux*(1-math.cos(ang)))+(uz*math.sin(ang))
   e=(math.cos(ang))+(uy*uy*(1-math.cos(ang)))
   f=(uy*uz*(1-math.cos(ang)))-(ux*math.sin(ang))
   g=(uz*ux*(1-math.cos(ang)))-(uy*math.sin(ang))
   h=(uz*uy*(1-math.cos(ang)))+(ux*math.sin(ang))
   i=math.cos(ang)+(uz*uz*(1-math.cos(ang)))
   bigmat=([[a,b,c],[d,e,f,],[g,h,i]])
   vect=np.dot(bigmat,vect3)
   return vect

def calcdist(a,b,carts):
   return np.linalg.norm(np.subtract(carts[a],carts[b]))

def elementID(massno):
   if massno < len(periodictable): return periodictable[massno]
   else: return "XX"

def bondiRadius(massno):
   #Bondi van der Waals radii for all atoms from: Bondi, A. J. Phys. Chem. 1964, 68, 441-452, except hydrogen, which is taken from Rowland, R. S.; Taylor, R. J. Phys. Chem. 1996, 100, 7384-7391
   #Radii that are not available in either of these publications have RvdW = 2.00 Angstrom

   bondi = [0.0,1.09, 1.40, 1.82,2.00,2.00,1.70,1.55,1.52,1.47,1.54,2.27,1.73,2.00,2.10,1.80,1.80,1.75,1.88,2.75,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.40,1.39,1.87,2.00,1.85,1.90,
            1.85,2.02,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.63,1.72,1.58,1.93,2.17,2.00,2.06,1.98,2.16,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.72,1.66,1.55,1.96,2.02,2.00,2.00,2.00,
            2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,2.00,1.86]
   if massno<len(bondi): radius = bondi[massno]
   else: radius = 2.0
   return radius


def calcopposite(atom1,atom2,angle,molcart):
   h=calcdist(atom1,atom2,molcart)
   d=h*math.sin(angle)
   return d
def calcadj(atom1,atom2,angle,molcart):
   h=calcdist(atom1,atom2,molcart)
   d=h*math.cos(angle)
   return d

def getcoords(atom,molcart):
   coords=[]
   for i in range(3):
      coords.append(molcart[atom][i])
   return coords

def avpoints(atomnos,molcart):
   xcoords=[]
   ycoords=[]
   zcoords=[]
   for a in atomnos:
      xcoords.append(molcart[a][0])
      ycoords.append(molcart[a][1])
      zcoords.append(molcart[a][2])
   syslength=len(xcoords)
   x=0;y=0;z=0
   for i in range(syslength):
      x=x+xcoords[i]
      y=y+ycoords[i]
      z=z+zcoords[i]
   x=x/syslength; y=y/syslength; z=z/syslength
   return round(x,8),round(y,8),round(z,8)

def distcalc(atom1,atom2):
   x=atom1[0]-atom2[0]
   y=atom1[1]-atom2[1]
   z=atom1[2]-atom2[2]
   dist = (x**2+y**2+z**2)**0.5
   return dist

def dprod(v1, v2): return sum((a*b) for a, b in zip(v1, v2))

def length(v): return math.sqrt(dprod(v, v))

def angle(v1, v2):
   val = dprod(v1, v2) / length(v1) / length(v2)
   if val > 0.999999: val = 1.0
   if val < -0.999999: val = -1.0
   return math.acos(val)

def dihedral(atoma,atomb,atomc,atomd):
   x1=atoma[0]
   y1=atoma[1]
   z1=atoma[2]
   x2=atomb[0]
   y2=atomb[1]
   z2=atomb[2]
   x3=atomc[0]
   y3=atomc[1]
   z3=atomc[2]
   x4=atomd[0]
   y4=atomd[1]
   z4=atomd[2]
   ax= (y2-y1)*(z2-z3)-(z2-z1)*(y2-y3)
   ay= (z2-z1)*(x2-x3)-(x2-x1)*(z2-z3)
   az= (x2-x1)*(y2-y3)-(y2-y1)*(x2-x3)
   bx= (y3-y2)*(z3-z4)-(z3-z2)*(y3-y4)
   by= (z3-z2)*(x3-x4)-(x3-x2)*(z3-z4)
   bz= (x3-x2)*(y3-y4)-(y3-y2)*(x3-x4)
   nbx= (y2-y3)*(z4-z3)-(z2-z3)*(y4-y3)
   nby= (z2-z3)*(x4-x3)-(x2-x3)*(z4-z3)
   nbz= (x2-x3)*(y4-y3)-(y2-y3)*(x4-x3)
   torsion=180.0/math.pi*math.acos((ax*bx+ay*by+az*bz)/(math.sqrt(ax*ax+ay*ay+az*az)*math.sqrt(bx*bx+by*by+bz*bz)))
   sign=180.0/math.pi*math.acos((nbx*(x2-x1)+nby*(y2-y1)+nbz*(z2-z1))/(math.sqrt(nbx*nbx+nby*nby+nbz*nbz)*math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1))))
   if sign<90.0:
      torsion=torsion*-1.0
   return torsion


#Read molecule data from an input file - currently Gaussian *com and *pdb supported
class getinData:
    def __init__(self, file):
        def getJOBTYPE(self, inlines):
            if fileformat == "com":
                for i in range(0,len(inlines)):
                    if inlines[i].find("#") > -1: self.JOBTYPE = inlines[i]
        
        def accepted_file(file):
            filesplit = file.split(".")
            if len(filesplit) > 1: # there is an extension
                for suffix in ["com", "pdb"]: # authorized suffix in ccparse
                    if filesplit[len(filesplit)-1] == suffix: return True
            return False
        
        def getCHARGE(self, inlines):
            if fileformat == "com":
                for i in range(0,len(inlines)):
                    if inlines[i].find("#") > -1:
                        if len(inlines[i+1].split()) == 0:
                            self.CHARGE = inlines[i+4].split()[0]
                            self.MULT = inlines[i+4].split()[1]
                        if len(inlines[i+2].split()) == 0:
                            self.CHARGE = inlines[i+5].split()[0]
                            self.MULT = inlines[i+5].split()[1]
        
        def getMEMREQ(self, inlines):
            if fileformat == "com":
                for i in range(0,len(inlines)):
                    if inlines[i].find("%mem") > -1: self.MEMREQ = inlines[i].split("=")[1].rstrip("\n")
        
        
        def getNPROC(self, inlines):
            if fileformat == "com":
                for i in range(0,len(inlines)):
                    if inlines[i].find("%nproc") > -1: self.NPROC = inlines[i].split("=")[1].rstrip("\n")
        
        
        def getATOMTYPES(self, inlines):
            if fileformat == "com":
                for i in range(0,len(inlines)):
                    if inlines[i].find("#") > -1:
                        if len(inlines[i+1].split()) == 0: start = i+5
                        if len(inlines[i+2].split()) == 0: start = i+6
                        break
                
                self.ATOMTYPES = []
                self.LEVELTYPES = []
                self.MMTYPES = []
                for i in range(start,len(inlines)):
                    if len(inlines[i].split()) ==0: break
                    else:
                        atominfo = inlines[i].split()[0]
                        atominfo = atominfo.split("-")[0]
                        if len(inlines[i].split()[0].split(atominfo))>1:
                            mminfo = inlines[i].split()[0].lstrip(atominfo)
                            self.MMTYPES.append(mminfo)
                        
                        self.ATOMTYPES.append(atominfo)
                        level = ""
                        for oniomlevel in ["H", "M", "L"]:
                            if inlines[i][4:].rfind(oniomlevel)>1:
                                level = inlines[i][4:][inlines[i][4:].rfind(oniomlevel):].rstrip("\n")
                        self.LEVELTYPES.append(level)
            if fileformat == "pdb":
                self.ATOMTYPES = []
                for i in range(0,len(inlines)):
                    if inlines[i].find("ATOM") > -1:
                       self.ATOMTYPES.append(int(inlines[i].split()[1]))
                    if inlines[i].find("HETATM")>-1:
                       self.ATOMTYPES.append(inlines[i].split()[-1])
        
        def getCONNECTIVITY(self, inlines, natoms):
            if fileformat == "com":
                for i in range(0,len(inlines)):
                    if inlines[i].find("#") > -1:
                        if len(inlines[i+1].split()) == 0: start = i+natoms+6
                        if len(inlines[i+2].split()) == 0: start = i+natoms+7
                        break
                
                if start < len(inlines):
                    self.CONNECTIVITY = []
                    j = 1
                    for i in range(start,len(inlines)):
                        if len(inlines[i].split()) != 0:
                            try: num = int(inlines[i].split()[0])
                            except ValueError: num = 0
                        if num == j:
                            bond=[]
                            neighbors=(len(inlines[i].split())-1)/2
                            if neighbors!=0:
                                for k in range(0,neighbors): bond.append((inlines[i].split()[1+2*k])+"__"+(inlines[i].split()[2+2*k]))
                            self.CONNECTIVITY.append(bond)
                            j = j+1
                    
                    if len(self.CONNECTIVITY) == natoms:
                        for i in range(0, natoms):
                            for partner in self.CONNECTIVITY[i]:
                                info = partner.split("__")
                                nextatom = int(info[0])-1
                                bondorder = float(info[1])
                                nope=0
                                for otherpartner in self.CONNECTIVITY[nextatom]:
                                    otherinfo = otherpartner.split("__")
                                    othernextatom = int(otherinfo[0])-1
                                    if othernextatom==i: nope=nope+1
                                if nope==0: self.CONNECTIVITY[nextatom].append(str(i+1)+"__"+info[1])
                    
                    self.OPTIONAL = []
                    for i in range(start+j,len(inlines)):
                        if len(inlines[i].split()) != 0: self.OPTIONAL.append(inlines[i])
            if fileformat == "pdb":
                self.CONNECTIVITY = []
                for i in range(self.NATOMS):
                    self.CONNECTIVITY.append([]) #pre-fill the list
                for i in range(0,len(inlines)):
                    if inlines[i].find("CONECT") > -1: # a connectivity
                        conect = inlines[i].split()
                        if len(conect) > 1: # contains at least a bond connectivity
                            try: num = int(conect[1])-1
                            except ValueError: num = 0
                            bond=[]
                            neighbors=(len(conect)-2)
                            if neighbors!=0:
                                for k in range(0,neighbors): bond.append((conect[2+k])+"__1.0") #1.0 default value for bond order: single bond
                            self.CONNECTIVITY[num].extend(bond) #add the new bonds
                # several bonds might have been added several time if the software was stupid. Lets clean
                for i in range(self.NATOMS):
                    self.CONNECTIVITY[i] = list(set(self.CONNECTIVITY[i]))
                #At the end, we should have a connectivity for each atom. We check the data is consistent and we fix it.
                if len(self.CONNECTIVITY) == natoms:
                    for i in range(0, natoms):
                        for partner in self.CONNECTIVITY[i]:
                            info = partner.split("__")
                            nextatom = int(info[0])-1
                            bondorder = float(info[1])
                            nope=0
                            for otherpartner in self.CONNECTIVITY[nextatom]:
                                otherinfo = otherpartner.split("__")
                                othernextatom = int(otherinfo[0])-1
                                if othernextatom==i: nope=nope+1
                            if nope==0: self.CONNECTIVITY[nextatom].append(str(i+1)+"__"+info[1]) # add the connectivity to the nextatom
                else:
                    print "Error: Number of connectivity inconsistent with coordinates"
                    
        def getCARTESIANS(self, inlines, natoms):
            if fileformat == "com":
                for i in range(0,len(inlines)):
                    if inlines[i].find("#") > -1:
                        if len(inlines[i+1].split()) == 0: start = i+5
                        if len(inlines[i+2].split()) == 0: start = i+6
                        break
                
                self.CARTESIANS = []
                for i in range(start,len(inlines)):
                    if len(inlines[i].split()) == 0: break
                    elif len(inlines[i].split()) == 4: self.CARTESIANS.append([float(inlines[i].split()[1]), float(inlines[i].split()[2]), float(inlines[i].split()[3])])
                    elif len(inlines[i].split()) > 4: self.CARTESIANS.append([float(inlines[i].split()[2]), float(inlines[i].split()[3]), float(inlines[i].split()[4])])
            if fileformat == "pdb":
                self.CARTESIANS = []
                for i in range(0,len(inlines)):
                    if inlines[i].find("ATOM") > -1:
                       self.CARTESIANS.append(float(inlines[i].split()[2:4]))
                    if inlines[i].find("HETATM")>-1:
                       self.CARTESIANS.append([float(inlines[i].split()[-6]), float(inlines[i].split()[-5]), float(inlines[i].split()[-4])])

        def getCONSTRAINED(self, optional):
            if fileformat == "com":
                self.CONSTRAINED = []
                for line in optional:
                    if line.find("X") > -1 and line.find("F") > -1: self.CONSTRAINED.append([int(line.split(" ")[1])-1])
                    if line.find("B") > -1 and line.find("F") > -1: self.CONSTRAINED.append([int(line.split(" ")[1])-1,int(line.split(" ")[2])-1])
                    if line.find("A") > -1 and line.find("F") > -1: self.CONSTRAINED.append([int(line.split(" ")[1])-1,int(line.split(" ")[2])-1]),int(line.split(" ")[3])-1
                    if line.find("D") > -1 and line.find("F") > -1: self.CONSTRAINED.append([int(line.split(" ")[1])-1,int(line.split(" ")[2])-1, int(line.split(" ")[3])-1, int(line.split(" ")[4])-1])
        
        if accepted_file(file):
            # default values
            self.CHARGE = 0
            self.MULT = 0
            self.MEMREQ = ""
            self.NPROC = ""
            self.ATOMTYPES = []
            self.LEVELTYPES = []
            self.MMTYPES = []
            self.CONNECTIVITY = []
            self.CARTESIANS = []
            self.CONSTRAINED = []
            # analyze
            filesplit = file.split(".")
            fileformat = filesplit[len(filesplit)-1]
            infile = open(file,"r")
            inlines = infile.readlines()
            self.NAME = file
            getJOBTYPE(self, inlines)
            getCHARGE(self, inlines)
            getMEMREQ(self, inlines)
            getNPROC(self, inlines)
            getATOMTYPES(self, inlines)
            self.NATOMS=len(self.ATOMTYPES)
            getCARTESIANS(self, inlines, self.NATOMS)
            getCONNECTIVITY(self, inlines, self.NATOMS)
            if hasattr(self, "OPTIONAL"): getCONSTRAINED(self, self.OPTIONAL)
            if len(self.ATOMTYPES) == 0 or len(self.CARTESIANS) ==0: print "\nFATAL ERROR: Input file [ %s ] cannot be read"%file
            #print "Input file data [ %s ]\n"% (file)
            #for i in range(len(self.CARTESIANS)):
            #    print "%3s  %3s     %8.3f %8.3f %8.3f      CONNECT %s\n"% (i+1, self.ATOMTYPES[i], self.CARTESIANS[i][0],self.CARTESIANS[i][1],self.CARTESIANS[i][2], self.CONNECTIVITY[i]) 
        else: print "\nError: Input file [ %s ] is not supported. [com, pdb] "

#Read molecule data from an output file #######################
class getoutData:
    def __init__(self, file):
        
        self.NAME = file
        
        def accepted_file(file):
            filesplit = file.split(".")
            if len(filesplit) > 1: # there is an extension
                for suffix in ["out", "log"]: # authorized suffix
                    if filesplit[len(filesplit)-1] == suffix: return True
            return False

        def getFORMAT(self, outlines):
            for i in range(0,len(outlines)):
                if outlines[i].find("MOPAC") > -1: self.FORMAT = "Mopac"; break
                if outlines[i].find("Gaussian") > -1: self.FORMAT = "Gaussian"; break
        
        def getCHARGE(self, outlines, format):
            if format == "Mopac":
                for i in range(0,len(outlines)):
                    if outlines[i].find("CHARGE ON SYSTEM") > -1:
                        self.CHARGE = int(outlines[i].split()[5])
                        self.MULT = 1
                    if outlines[i].find("STATE CALCULATION") > -1:
                        if outlines[i].split()[0] == "SINGLET": self.MULT = 1
                        if outlines[i].split()[0] == "DOUBLET": self.MULT = 2
                        if outlines[i].split()[0] == "TRIPLET": self.MULT = 3
                        break
            if format == "Gaussian":
                for i in range(0,len(outlines)):
                    if outlines[i].find("Charge = ") > -1:
                        self.CHARGE = int(outlines[i].split()[2])
                        self.MULT = int(outlines[i].split()[5].rstrip("\n"))
                        break
        
        
        def getATOMTYPES(self, outlines, format):
            self.ATOMTYPES = []
            self.CARTESIANS = []
            if format == "Mopac":
                for i in range(0,len(outlines)):
                    if outlines[i].find("CHEMICAL") > -1: standor = i+3
                    if outlines[i].find("Empirical Formula") > -1:
                        self.NATOMS = int((outlines[i].split("="))[1].split()[0])
                
                if hasattr(self, "NATOMS"):
                    for i in range (standor,standor+self.NATOMS):
                        outlines[i] = outlines[i].replace("*", " ")
                        self.ATOMTYPES.append((outlines[i].split()[3]))
                        self.CARTESIANS.append([float(outlines[i].split()[5]), float(outlines[i].split()[6]), float(outlines[i].split()[7])])
            
            if format == "Gaussian":
                for i in range(0,len(outlines)):
                    if outlines[i].find("Standard orientation") > -1:
                        standor = i
                    if outlines[i].find("Rotational constants") > -1 and outlines[i-1].find("-------") > -1:
                        self.NATOMS = i-standor-6
                try: standor
                except NameError: pass
                else:
                    for i in range (standor+5,standor+5+self.NATOMS):
                        self.ATOMTYPES.append(elementID(int(outlines[i].split()[1])))
                        self.CARTESIANS.append([float(outlines[i].split()[3]),float(outlines[i].split()[4]),float(outlines[i].split()[5])])
        
        
        def getFREQS(self, outlines, format):
            self.FREQS = []
            if format == "Gaussian":
                for i in range(0,len(outlines)):
                    if outlines[i].find("Frequencies") > -1:
                        self.FREQS.append(float(outlines[i].split()[2]))
                        if len(outlines[i].split()) > 3: self.FREQS.append(float(outlines[i].split()[3]))
                        if len(outlines[i].split()) > 4: self.FREQS.append(float(outlines[i].split()[4]))
                if len(self.FREQS) > 0:
                    for i in range(0,len(outlines)):
                        if outlines[i].find("Zero-point correction") > -1: self.ZPE = float(outlines[i].split()[2])
                        if outlines[i].find("thermal Enthalpies") > -1: self.ENTHALPY = float(outlines[i].split()[6])
                        if outlines[i].find("thermal Free Energies") > -1: self.GIBBS = float(outlines[i].split()[7])
        
        def getCPU(self, outlines, format):
            days = 0; hours = 0; mins = 0; secs = 0
            if format == "Mopac":
                for i in range(0,len(outlines)):
                    if outlines[i].find("TOTAL CPU TIME") > -1:
                        secs = secs + int(float(outlines[i].split()[3]))
            if format == "Gaussian":
                for i in range(0,len(outlines)):
                    if outlines[i].find("Job cpu time") > -1:
                        days = days + int(outlines[i].split()[3])
                        hours = hours + int(outlines[i].split()[5])
                        mins = mins + int(outlines[i].split()[7])
                        secs = secs + int(float(outlines[i].split()[9]))
            self.CPU=[days,hours,mins,secs]
        
        
        
        def getENERGY(self, outlines, format):
            if format == "Mopac":
                for i in range(0,len(outlines)):
                    if outlines[i].find("TOTAL ENERGY") > -1:
                        self.ENERGY = 0.036749309*(float(outlines[i].split()[3])) # eV to hartree
            if format == "Gaussian":
                uff = 0
                am1 = 0
                pm3 = 0
                scf = 0
                oniom = 0
                for i in range(0,len(outlines)):
                    if outlines[i].find(" UFF") > -1: uff = i
                    if outlines[i] .find("AM1") > -1: am1 = i
                    if outlines[i].find("PM3") > -1: pm3 = i
                    if outlines[i].find("ONIOM") > -1: oniom = i
                    if outlines[i].find("SCF Done") > -1: scf = i
                
                calctype = [uff,am1,pm3,oniom,scf]
                for i in range(0,len(outlines)):
                    if scf == max(calctype) and outlines[i].find("SCF Done") > -1 and outlines[i].find("Initial convergence to 1.0D-05 achieved")==-1: # Get energy from HF or DFT calculation
                        self.ENERGY = (float(outlines[i].split()[4]))
                    if oniom == max(calctype) and outlines[i].find("ONIOM: extrapolated energy") > -1: # Get energy from ONIOM calculation
                        self.ENERGY = (float(outlines[i].split()[4]))
                    if pm3 == max(calctype) or am1 == max(calctype) or uff == max(calctype):
                        if outlines[i].find("Energy= ") > -1 and outlines[i].find("Predicted")==-1 and outlines[i].find("Thermal")==-1: # Get energy from Semi-empirical or Molecular Mechanics calculation
                            self.ENERGY = (float(outlines[i].split()[1]))
                    if outlines[i].find("Total free energy in solution") > -1:
                        self.SOLVENERGY = (float(outlines[i+1].split()[7]))

        
        if not os.path.exists(self.NAME):
            print ("\nFATAL ERROR: Output file [ %s ] does not exist"%file)
        else:
            if accepted_file(file): #check extensions
                outfile = open(self.NAME,"r")
                outlines = outfile.readlines()
                getFORMAT(self, outlines)
                if hasattr(self, "FORMAT"):
                    getCHARGE(self, outlines, self.FORMAT)
                    getENERGY(self, outlines, self.FORMAT)
                    getATOMTYPES(self, outlines, self.FORMAT)
                    getFREQS(self, outlines, self.FORMAT)
                    getCPU(self, outlines, self.FORMAT)
            else:
                print ("\nFATAL ERROR: Output file [ %s ] is not supported [.out, .log]"%file)
###############################################################

def concheck(conpar,val):
   cons=[]
   for a in range(len(conpar)):
      for b in range(len(conpar[a])):
         if val ==conpar[a][0]:
            for c in range(len(conpar[a])-1):
               cons.append(conpar[a][c+1])
            return cons

def twod_dist(a,b,c):
   vect1=np.subtract(a,b)
   vect2=np.subtract(b,c)
   ang=angle(vect1,vect2)
   return math.sin(ang)*np.linalg.norm(vect1)

def twod_vect(a,b,c):
   vect1=np.subtract(a,b)
   vect2=np.subtract(b,c)
   ang=angle(vect1,vect2)
   nvect2=vect2/np.linalg.norm(vect2)
   return ((math.cos(ang)*np.linalg.norm(vect1))*nvect2)+b

def twod_rot(vect,theta):
   a=math.cos(theta)
   b=math.sin(theta)
   mat=[[a,-b],[b,a]]
   vect=np.dot(mat,vect)
   return vect

# Generate Sterimol atom type from connectivity data
def generate_atom_types(atomtype, cn):
   st_types = []
   for i in range(0,len(atomtype)):
      atom = atomtype[i]
      if atom == "H": st_types.append("H")
      elif atom == "P": st_types.append("P")
      elif atom == "F": st_types.append("F")
      elif atom == "Cl": st_types.append("C1")
      elif atom == "Br": st_types.append("B1")
      elif atom == "I": st_types.append("I")
      elif atom == "O": #Sterimol distinguishes between "normal", and double-bonded O atoms
         if cn[i] < 1.5: st_types.append("O2")
         if cn[i] > 1.5: st_types.append("O")
      elif atom == "S": #Sterimol distinguishes between "normal", tetrahedral, and octohedral S atoms
         if cn[i] < 2.5: st_types.append("S")
         if 5.5 > cn[i] > 2.5: st_types.append("S4")
         if cn[i] > 5.5: st_types.append("S1")
      elif atom == "N": #Sterimol distinguishes between tetrahedral and planar (amide) N atoms
         if cn[i] > 2.5: st_types.append("N")
         if cn[i] < 2.5: st_types.append("C6/N6")
      elif atom == "C": #Sterimol distinguishes between myriad types of C atoms ...
         if cn[i] < 2.5: st_types.append("C3")
         if 3.5 > cn[i] > 2.5: # need to differentiate between sp2 carbon and aromatic carbon ...
            st_types.append("C6/N6") # assumes aromatic rather than sp2
         if cn[i] > 3.5: st_types.append("C")
   return st_types


# Calculation of atomic coordination numbers (taken from Grimme's DFTD3 definitions)
def ncoord(natom, rcov, atomtype, coords):
   max_elem = 94
   k1 = 16.0
   k2 = 4.0/3.0
   cn =[]
   for i in range(0,natom):
      xn = 0.0
      for iat in range(0,natom):
         if iat != i:
            dx = coords[iat][0] - coords[i][0]
            dy = coords[iat][1] - coords[i][1]
            dz = coords[iat][2] - coords[i][2]
            r2 = dx*dx+dy*dy+dz*dz
            r = math.pow(r2,0.5)
            r = r
            for k in range(0,max_elem):
               if atomtype[i].find(elements[k])>-1:Zi=k
               if atomtype[iat].find(elements[k])>-1:Ziat=k

            rco = rcov[Zi]+rcov[Ziat]
            rco = rco*k2
            rr=rco/r
            damp=1.0/(1.0+math.exp(-k1*(rr-1.0)))
            xn=xn+damp
      cn.append(xn)
   return cn

def linearcheck(carts):
   ans=0;xgrad=[];ygrad=[]
   for row in carts:xgrad.append(round(np.gradient(row)[0],4));ygrad.append(round(np.gradient(row)[1],4))
   if min(xgrad)==max(xgrad) and min(ygrad)==max(ygrad):ans=1
   return ans

class calcSterimol:
   def __init__(self, file, radii, atomA, atomB,verbose):

      if len(file.split(".com"))>1 or len(file.split(".gjf"))>1: fileData = getinData(file)
      if len(file.split(".out"))>1 or len(file.split(".log"))>1: fileData = getoutData(file)

      # initialize the array of atomic vdw radii
      molcart = fileData.CARTESIANS; atomtype = fileData.ATOMTYPES; natoms = len(molcart); vdw_radii = []

      if radii == "cpk":
         atomic_co_no = ncoord(natoms, rcov, atomtype, molcart)
         sterimol_types = generate_atom_types(atomtype, atomic_co_no)
         #print sterimol_types
         for i in range(0,natoms):
            for j in range(0,len(sterimol_atomtypes)):
               if sterimol_types[i] == sterimol_atomtypes[j]: vdw_radii.append(cpk_radii[j]/100.00)

      if radii == "bondi":
         for i in range(0,natoms): vdw_radii.append(bondiRadius(periodictable.index(fileData.ATOMTYPES[i])))

# Define vector along the L-axis connecting base atom and the next attached atom
# subtract one since the array starts from zero not one
      atomA = atomA - 1; atomB = atomB - 1
      next_atom = molcart[atomB]
      vect1=np.subtract(getcoords(atomA,molcart),next_atom)
      if verbose == True:
          print "   Atoms", atomA, "and", atomB, "define the L-axis and direction", vect1

          print "\n", "   Atom ".ljust(9), "  Xco/A".rjust(9), "  Yco/A".rjust(9), "  Zco/A".rjust(9), " VdW/pm".rjust(9)
          print "   ##############################################"
      # Remove the base atom from the list of atoms to be considered for sterics (after printing all)
      atomlist = list(xrange(0,natoms))
      if verbose == True:
          for atom in atomlist:
             if radii == "cpk": print "  ", sterimol_types[atom].ljust(6),
             if radii == "bondi": print "  ", atomtype[atom].ljust(6),
             for coord in molcart[atom]:
                if coord < 0.0: print "   %.3f".rjust(6) % coord,
                else: print "    %.3f".rjust(6) % coord,
             print "    %.1f" % round(vdw_radii[atom]*100)
      atomlist.remove(atomA)

      adjlist=[]; opplist=[]; theta=[]
      for i in atomlist:
         vect2=np.subtract(getcoords(atomA,molcart),getcoords(i,molcart))
         oppdist=calcopposite(atomA,i,angle(vect1,vect2),molcart)
         opplist.append(oppdist+vdw_radii[i])
         adjdist=calcadj(atomA,i,angle(vect1,vect2),molcart)
         #minadjlist.append(adjdist-vdw_radii[i])
         adjlist.append(adjdist+vdw_radii[i])

      B5=max(opplist)
   #self.lval=max(adjlist)-minval
   # A bit weird, but seems like original sterimol adds on the difference between the bond length and vdw radius of atom B. For a C-H bond this is 1.50 - 1.10 = 0.40 Angstrom)
      self.lval=max(adjlist)+0.40

      ###Useful - do not delete!
      #print "   B5 atom", atomlist[opplist.index(max(opplist))]+1, "distance", max(opplist)
      #print "   Highest atom", atomlist[adjlist.index(max(adjlist))]+1,"distance", max(adjlist),"\n   Lowest atom", atomlist[minadjlist.index(min(minadjlist))]+1,"distance", min(minadjlist)

      zcarts=[]#zeroed carts
      for i in atomlist: zcarts.append(np.subtract(molcart[i],molcart[atomA]))
      zvect=[0,0,1]
      zcent=np.subtract(next_atom,molcart[atomA])
      for cart in range(len(zcarts)):
         zcoord= rotrel(zcent,zvect,zcarts[cart])
         zcarts[cart]=zcoord
      twodcarts=[]
      for row in zcarts: twodcarts.append([row[0],row[1]])
      fragrad=[]#radii of fragment atoms
      for t in atomlist: fragrad.append(vdw_radii[t])
      singledist=[]
      for t in range(len(fragrad)):
         d=np.linalg.norm(twodcarts[t])#;print d
         d=d+fragrad[t]
         singledist.append(d)
      self.newB5=max(singledist) #This is the same as the 3D calculated value from above

      center=[0,0]
      vlist=[]#list of distances from the origin to the tangential vectors
      alist=[]#list of atoms between which the tangential vectors pass through no other atoms
      iav=[]#interatomic vectors
      sym=symcheck(twodcarts)
      for x in range(len(twodcarts)):
         if sym==1:
            twodcarts[x][0]=twodcarts[x][0]+0.000001
            twodcarts[x][1]=twodcarts[x][1]+0.000001
         for y in range(len(twodcarts)):
            if x!=y:
               try:nvect= (twod_vect(center,twodcarts[x],twodcarts[y]))#origin normal vector to connecting atomic centers vector
               except ValueError:nvect=[0,0]
               iav=np.subtract(twodcarts[x],twodcarts[y])#interatomic vector
               iad=np.linalg.norm(iav)#interatomic distance
               try:theta=math.asin((fragrad[y]-fragrad[x])/iad)#calculates angle by which to rotate vdw radii before adding
               except ValueError: theta=np.pi/2
               try:unvect=nvect/np.linalg.norm(nvect)
               except RuntimeWarning:pass#unvect=[0,0]
               xradv=twod_rot(unvect*fragrad[x],theta)
               yradv=twod_rot(unvect*fragrad[y],theta)
               mvect= (twod_vect(center,twodcarts[x]-xradv,twodcarts[y]-yradv))
               nvect= (twod_vect(center,twodcarts[x]+xradv,twodcarts[y]+yradv))#origin normal vector to connecting atomic surfaces tangential vector
               newx=twodcarts[x]+xradv
               newy=twodcarts[y]+yradv
               mewx=twodcarts[x]-xradv
               mewy=twodcarts[y]-yradv
               if np.cross(nvect,xradv)<0.000000001 and theta!=np.pi/2:
                  satpoint=[]#Satisfied points not within range of tangential vector
                  for z in range(len(twodcarts)):
                     pvdist=twod_dist(twodcarts[z],newx,newy)
                     if z!=x and z!=y and pvdist>(fragrad[z]-0.0001):satpoint.append(pvdist)
                  if len(satpoint)==len(atomlist)-2:vlist.append(np.linalg.norm(nvect));alist.append([x,y]);#print x,y
                  satpoint=[]
                  for z in range(len(twodcarts)):
                     pvdist=twod_dist(twodcarts[z],mewx,mewy)
                     if z!=x and z!=y and pvdist>(fragrad[z]-0.0001):satpoint.append(pvdist)
                  if len(satpoint)==len(atomlist)-2:vlist.append(np.linalg.norm(mvect));alist.append([x,y])
      if linearcheck(twodcarts)==1:self.B1 = max(fragrad)
      elif len(vlist) > 0: self.B1=min(vlist)
      else: self.B1 = max(fragrad)

def symcheck(carts):#Add symmetry criteria
   center=[0,0]
   distlist=[]
   distlist.append(10)
   for a in range(len(carts)):
      for b in range(len(carts)):
         if a!=b:
            dist=np.linalg.norm(twod_vect(center,carts[a],carts[b]))
            distlist.append(dist)
   if min(distlist)<0.0000000001:ans=1
   else:ans=0
   return ans

def calcSandwich(file):
   metalatoms=[]
   if file.split(".")[1]=="log" or file.split(".")[1]=="out":fileData=getoutData(file.split(".")[0])
   if file.split(".")[1]=="com" or file.split(".")[1]=="gjf":fileData=getinData(file.split(".")[0])
   for i in range(len(fileData.ATOMTYPES)):
      if fileData.ATOMTYPES[i] in metals:metalatoms.append(i)

   ivals=[]
   jvals=[]
   for i in range(len(fileData.ATOMTYPES)):
      for j in range(len(fileData.ATOMTYPES)):
         dist = ((fileData.CARTESIANS[i][0]-fileData.CARTESIANS[j][0])**2 +(fileData.CARTESIANS[i][1]-fileData.CARTESIANS[j][1])**2+(fileData.CARTESIANS[i][2]-fileData.CARTESIANS[j][2])**2)**0.5
         if 0.01<dist<1.511 and fileData.ATOMTYPES[j] == "C" and fileData.ATOMTYPES[i] == "C":
            ivals.append(i)
            jvals.append(j)
   conpar=[]
   for a in range(len(ivals)):
      rar=[]
      rar.append(ivals[a])

      for b in range(len(ivals)):
         if ivals[a]==ivals[b]:rar.append(jvals[b])
      if rar not in conpar:conpar.append(rar)

   allrings=[]
   for a in range(len(conpar)):
      z=conpar[a][0]
      for b in concheck(conpar,z):
         y=b
         for c in concheck(conpar,y):
            x=c
            for d in concheck(conpar,x):
               w=d
               for e in concheck(conpar,w):
                  v=e
                  rar=[]
                  rar.extend([z,y,x,w,v])
                  if z in concheck(conpar,v) and sorted(rar) not in allrings and len(set(rar))==5:allrings.append(sorted(rar))
                  for f in concheck(conpar,v):
                     u=f
                     tar=[]
                     tar.extend([z,y,x,w,v,u])
                     if z in concheck(conpar,u) and sorted(tar) not in allrings and len(set(tar))==6:allrings.append(sorted(tar))



   if not allrings:
      for ma in metalatoms:
         for s in range(len(fileData.CARTESIANS)):
            if 0.1<np.linalg.norm(np.subtract(fileData.CARTESIANS[ma],fileData.CARTESIANS[s]))<2.1:allrings.append([s,s,s,s,s])

   mcdists=[]
   mcdist=9999
   for ring in allrings:

      if len(ring)==5:
         tolman=[]
         cent=avpoints(ring,fileData.CARTESIANS)
         m=fileData.CARTESIANS[metalatoms[0]]
         tempmcdist=mcdist
         mcdist=distcalc(m,cent)
         for b in metalatoms:#find closest metal to ring
            m=fileData.CARTESIANS[b]
            if mcdist>=distcalc(m,cent):mcdist=distcalc(m,cent);metal=b
         mcdists.append([mcdist,metal])
         frag=getfragment(ring[0],fileData.CARTESIANS)
         vect1=np.subtract(getcoords(metal,fileData.CARTESIANS),cent)
         if tempmcdist==mcdist:break#Stops if dealing with identical ring system as before (intended for symmetric dimers)
         adjlist=[]
         minadjlist=[]
         opplist=[]
         alpha=[]
         beta=[]
         theta=[]#Candidate Tolman angle substituent
         omega=[]#standardised atom "dihedral" orientation
         ringang=[]
         for i in frag:
            vect2=np.subtract(getcoords(metal,fileData.CARTESIANS),getcoords(i,fileData.CARTESIANS))
            oppdist=calcopposite(metal,i,angle(vect1,vect2),fileData.CARTESIANS)
            opplist.append(oppdist+genradii(i,fileData.CARTESIANS,fileData.ATOMTYPES))
            adjdist=calcadj(metal,i,angle(vect1,vect2),fileData.CARTESIANS)
            minadjlist.append(adjdist-genradii(i,fileData.CARTESIANS,fileData.ATOMTYPES))
            adjlist.append(adjdist+genradii(i,fileData.CARTESIANS,fileData.ATOMTYPES))
            alpha.append(angle(vect1,vect2))
            hyp=distcalc(getcoords(i,fileData.CARTESIANS),getcoords(metal,fileData.CARTESIANS))
            beta.append(math.asin(genradii(i,fileData.CARTESIANS,fileData.ATOMTYPES)/hyp))
            theta.append(alpha[-1]+beta[-1])
            if ring[0]!=ring[1]:omega.append(dihedral([10,10,10],getcoords(metal,fileData.CARTESIANS),cent,getcoords(i,fileData.CARTESIANS)))
            if ring[0]!=ring[1] and i in ring:ringang.append(dihedral([10,10,10],getcoords(metal,fileData.CARTESIANS),cent,getcoords(i,fileData.CARTESIANS)))
         B5=max(opplist)#Bondi
         lval=max(adjlist)-min(minadjlist)
         interval=180/len(ring)
         if ring[0]!=ring[1]:
            for k in ringang:
               tlist=[];tang=[];tfrag=[]

               for h in range(len(frag)):
                  if k-interval<omega[h]<k+interval:tlist.append(frag[h])
                  if k>(180-interval) and k-interval<omega[h]+360<k+interval:tlist.append(frag[h])
                  if k<-(180-interval) and k-interval<omega[h]-360<k+interval:tlist.append(frag[h])
               for t in range(len(frag)):
                  if frag[t] in tlist: tang.append(theta[t]);tfrag.append(frag[t])
               tolman.append(math.degrees(max(tang)))
            x=0
            for c in tolman:
               x=x+c
            tolmanCA=round(2*(x/len(tolman)),3)
         else:tolmanCA=0
         smcdist=round(mcdist,3);lval=round(lval,3);lval=round(lval,3);B5=round(B5,3)
         molcart=fileData.CARTESIANS
         zcarts=[]
         for i in frag:
            zcarts.append(np.subtract(molcart[i],molcart[metal]))
         zvect=[0,0,1]
         zcent=np.subtract(cent,molcart[metal])
         for cart in range(len(zcarts)):
            zcoord= rotrel(zcent,zvect,zcarts[cart])
            zcarts[cart]=zcoord
         twodcarts=[]
         for row in zcarts:
            twodcarts.append([row[0],row[1]])
         fragrad=[]#radii of fragment atoms
         for t in frag:
            fragrad.append(genradii(t,fileData.CARTESIANS,fileData.ATOMTYPES))
         singledist=[]
         for t in range(len(fragrad)):
            d=np.linalg.norm(twodcarts[t])#;print d
            d=d+fragrad[t]
            singledist.append(d)
         newB5=round(max(singledist),3)#This is the same as the 3D calculated value from above

         center=[0,0]
         vlist=[]#list of distances from the origin to the tangential vectors
         alist=[]#list of atoms between which the tangential vectors pass through no other atoms
         iav=[]#interatomic vectors

         for x in range(len(twodcarts)):
            for y in range(len(twodcarts)):
               if x!=y:
                  try:nvect= (twod_vect(center,twodcarts[x],twodcarts[y]))#origin normal vector to connecting atomic centers vector
                  except ValueError:nvect=[0,0]
                  iav=np.subtract(twodcarts[x],twodcarts[y])#interatomic vector
                  iad=np.linalg.norm(iav)#interatomic distance
                  try:theta=math.asin((fragrad[y]-fragrad[x])/iad)#calculates angle by which to rotate vdw radii before adding
                  except ValueError: theta=np.pi/2
                  try:unvect=nvect/np.linalg.norm(nvect)
                  except RuntimeWarning:pass#unvect=[0,0]
                  xradv=twod_rot(unvect*fragrad[x],theta)
                  yradv=twod_rot(unvect*fragrad[y],theta)
                  nvect= (twod_vect(center,twodcarts[x]+xradv,twodcarts[y]+yradv))#origin normal vector to connecting atomic surfaces tangential vector
                  newx=twodcarts[x]+xradv
                  newy=twodcarts[y]+yradv
                  if np.cross(nvect,xradv)<0.000000001 and theta!=np.pi/2:
                     satpoint=[]#Satisfied points not within range of tangential vector
                     for z in range(len(twodcarts)):
                        pvdist=twod_dist(twodcarts[z],newx,newy)
                        if z!=x and z!=y and pvdist>fragrad[z]:satpoint.append(pvdist)
                     if len(satpoint)==len(frag)-2:vlist.append(np.linalg.norm(nvect));alist.append([x,y])#;print x,y
         B1=round(min(vlist),3)
         print "   "+file.ljust(25),str(tolmanCA).rjust(9), str(smcdist).rjust(9), str(lval).rjust(9),str(B1).rjust(9), str(newB5).rjust(9)

molmod=[['Bq', 0, 0, 0, 0],
        ['H', 1, 1, 1, 1],
        ['He', 0, 0, 0, 0],
        ['Li', 0, 0, 0, 0],
        ['Be', 0, 0, 0, 0],
        ['B', 0, 0, 0, 0],
        ['C', 0, 1.6, 1.6, 1.5],
        ['N', 1.45, 1.45, 1.5, 1.25],
        ['O', 1.35, 1.35, 1.35, 0],
        ['F', 1.35, 1.35, 0, 0],
        ['Ne', 0, 0, 0, 0],
        ['Na', 0, 0, 0, 0],
        ['Mg', 0, 0, 0, 0],
        ['Al', 0, 0, 0, 0],
        ['Si', 2.1, 2.1, 2.1, 2.1],
        ['P', 0, 0, 0, 0],
        ['S', 0, 0, 0, 0],
        ['Cl', 1.8, 0, 0, 0],
        ['Ar', 0, 0, 0, 0],
        ['K', 0, 0, 0, 0],
        ['Ca', 0, 0, 0, 0],
        ['Sc', 0, 0, 0, 0],
        ['Ti', 0, 0, 0, 0],
        ['V', 0, 0, 0, 0],
        ['Cr', 0, 0, 0, 0],
        ['Mn', 0, 0, 0, 0],
        ['Fe', 0, 0, 0, 0],
        ['Co', 0, 0, 0, 0],
        ['Ni', 0, 0, 0, 0],
        ['Cu', 0, 0, 0, 0],
        ['Zn', 0, 0, 0, 0],
        ['Ga', 0, 0, 0, 0],
        ['Ge', 0, 0, 0, 0],
        ['As', 0, 0, 0, 0],
        ['Se', 0, 0, 0, 0],
        ['Br', 1.95, 0, 0, 0],
        ['Kr', 0, 0, 0, 0],
        ['Rb', 0, 0, 0, 0],
        ['Sr', 0, 0, 0, 0],
        ['Y', 0, 0, 0, 0],
        ['Zr', 0, 0, 0, 0],
        ['Nb', 0, 0, 0, 0],
        ['Mo', 0, 0, 0, 0],
        ['Tc', 0, 0, 0, 0],
        ['Ru', 0, 0, 0, 0],
        ['Rh', 0, 0, 0, 0],
        ['Pd', 0, 0, 0, 0],
        ['Ag', 0, 0, 0, 0],
        ['Cd', 0, 0, 0, 0],
        ['In', 0, 0, 0, 0],
        ['Sn', 0, 0, 0, 0],
        ['Sb', 0, 0, 0, 0],
        ['Te', 0, 0, 0, 0],
        ['I', 2.15, 0, 0, 0],
        ['Xe', 0, 0, 0, 0],
        ['Cs', 0, 0, 0, 0],
        ['Ba', 0, 0, 0, 0],
        ['La', 0, 0, 0, 0],
        ['Ce', 0, 0, 0, 0],
        ['Pr', 0, 0, 0, 0],
        ['Nd', 0, 0, 0, 0],
        ['Pm', 0, 0, 0, 0],
        ['Sm', 0, 0, 0, 0],
        ['Eu', 0, 0, 0, 0],
        ['Gd', 0, 0, 0, 0],
        ['Tb', 0, 0, 0, 0],
        ['Dy', 0, 0, 0, 0],
        ['Ho', 0, 0, 0, 0],
        ['Er', 0, 0, 0, 0],
        ['Tm', 0, 0, 0, 0],
        ['Yb', 0, 0, 0, 0],
        ['Lu', 0, 0, 0, 0],
        ['Hf', 0, 0, 0, 0],
        ['Ta', 0, 0, 0, 0],
        ['W', 0, 0, 0, 0],
        ['Re', 0, 0, 0, 0],
        ['Os', 0, 0, 0, 0],
        ['Ir', 0, 0, 0, 0],
        ['Pt', 0, 0, 0, 0],
        ['Au', 0, 0, 0, 0],
        ['Hg', 0, 0, 0, 0],
        ['Tl', 0, 0, 0, 0],
        ['Pb', 0, 0, 0, 0],
        ['Bi', 0, 0, 0, 0],
        ['Po', 0, 0, 0, 0],
        ['At', 0, 0, 0, 0],
        ['Rn', 0, 0, 0, 0],
        ['Fr', 0, 0, 0, 0],
        ['Ra', 0, 0, 0, 0],
        ['Ac', 0, 0, 0, 0],
        ['Th', 0, 0, 0, 0],
        ['Pa', 0, 0, 0, 0],
        ['U', 0, 0, 0, 0],
        ['Np', 0, 0, 0, 0],
        ['Pu', 0, 0, 0, 0],
        ['Am', 0, 0, 0, 0],
        ['Cm', 0, 0, 0, 0],
        ['Bk', 0, 0, 0, 0],
        ['Cf', 0, 0, 0, 0],
        ['Es', 0, 0, 0, 0],
        ['Fm', 0, 0, 0, 0],
        ['Md', 0, 0, 0, 0],
        ['No', 0, 0, 0, 0],
        ['Lr', 0, 0, 0, 0],
        ['Rf', 0, 0, 0, 0],
        ['Db', 0, 0, 0, 0],
        ['Sg', 0, 0, 0, 0],
        ['Bh', 0, 0, 0, 0],
        ['Hs', 0, 0, 0, 0],
        ['Mt', 0, 0, 0, 0],
        ['Ds', 0, 0, 0, 0],
        ['Rg', 0, 0, 0, 0],
        ['Uub', 0, 0, 0, 0],
        ['Uut', 0, 0, 0, 0],
        ['Uuq', 0, 0, 0, 0],
        ['Uup', 0, 0, 0, 0],
        ['Uuh', 0, 0, 0, 0],
        ['Uus', 0, 0, 0, 0],
        ['Uuo', 0, 0, 0, 0],]

elements = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si",
            "P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni",
            "Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo",
            "Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba",
            "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
            "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
            "Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
            "Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq",
            "Uup","Uuh","Uus","Uuo"]

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
