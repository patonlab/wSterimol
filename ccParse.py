#!/usr/bin/python

###############################################################
#                         ccParse.py                          #
#                Reads compchem job file(s)                   #
###############################################################

#Python Libraries 
import subprocess, sys, os

## Check for integer when parsing ##
def is_number(s):
   try: int(s); return True
   except ValueError: return False

#Some useful arrays for chemists
periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr", "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl", "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo","Bq"]

def elementID(massno):
   if massno < len(periodictable): return periodictable[massno]
   else: return "XX"

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
               # print "%3s  %3s     %8.3f %8.3f %8.3f      CONNECT %s\n"% (i+1, self.ATOMTYPES[i], self.CARTESIANS[i][0],self.CARTESIANS[i][1],self.CARTESIANS[i][2], self.CONNECTIVITY[i]) 
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

