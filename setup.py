#!/usr/bin/python
from __future__ import print_function, absolute_import

import sys, os

class Setup:
    def __init__(self, log, path = "default", exe = "default"):
        #default values
        self.software = "MOPAC" # Mopac by default
        self.loaded = True # good by default
        self.exe = "C:\PROGRA~1\MOPAC\MOPAC2016.exe" #default install for Mopac in Windows
        self.SE = "PM6-DH2" # default force field
        self.charge = "0" # uncharged
        self.scf = "" # optimisation
        self.rmsd_cutoff = 0.1 # Angstroms
        self.memories = "8"
        self.procsshared = "8"
        self.leveloftheory = "wb97xd/6-31g(d)"
        self.spin = "1"
        self.RJCT = 0.50
        self.singlepointcalculation = ""
        self.Temperature = 298
        self.energywindow_cutoff = []
        self.print_cutoff = 5.0
        self.angle_count = 5
        self.radii = "cpk"
        if path != "default":
            self.path = os.path.join(path, "setup.ini")
        else:
            self.path = "setup.ini"
        #load setup.ini
        if os.path.exists(self.path):
            file = open(self.path,"r")
            lines = file.readlines()
            file.close()
            for k in range(len(lines)):
                config = lines[k].split('#')
                if len(config) > 1: # there is a commented line, just keep the first bit
                    lines[k] = config[0]
                config = lines[k].split('=')
                if len(config) > 2: # there is an "=" in the keyword
                    str = config[1]
                    for i in range(2,len(config)):
                        str = str + "=" + config[i]
                    config2 = []
                    config2.append(config[0])
                    config2.append(str)
                    config = config2
                if len(config) == 2: #if not, we ignore. It is an empty config value or a commentline
                    for i in range(len(config)):
                        config[i] = config[i].strip(' ').strip('\n').strip('\t')
                    if config[0] == "PROG_EXEC":
                        if exe != "default": #define a new value for exe
                            if os.path.exists(exe):
                                lines[k] = "PROG_EXEC = \'"+exe+"\'\n"
                                replacethisfile = open(self.path,"w")
                                for line in lines: replacethisfile.write(line)
                                replacethisfile.close()
                                config[1] = exe
                            else:
                                log.write("Warning: Specified path to executable doesn't exist. Use setup.ini value.")
                        if not os.path.exists(config[1]):
                            self.loaded = False
                            log.write("Error: Specified path to executable in self.ini doesn't exist. [%s]" % config[1])
                            return
                        self.exe = config[1]
                    elif config[0] == "SEMI_EMPIRICAL":
                        self.SE = config[1]
                    elif config[0] == "CHARGE":
                        if self.number(config[1]):
                            self.charge = int(config[1])
                        else: log.write("Warning: the charge is not a number [%s]. Default Neutral." % config[1])
                    elif config[0] == "OPTIMISATION": # single point calculation # yes (SCF) or no (optimisation)
                        config[1] = config[1].lower()
                        if config[1] == "no":
                            self.scf = "1SCF"
                        else:
                            self.scf = "" # empty is optimisation by default
                    elif config[0] == "RMSD_CLUSTER_OPT": # RMSD clusterisation cut-off for optimised structures
                        if self.number(config[1]) == True:
                            self.rmsd_cutoff = float(config[1])
                        else:
                            log.write("Warning: RMSD_CLUSTER_OPT value is not a number [%s]. Default value 0.1 Angstroms. %s " % (config[1], self.number(config[1])))
                    elif config[0] == "MEMORIES": # Memories
                        if self.number(config[1]) == True:
                            self.memories = int(config[1])
                        else:
                            log.write("Warning: MEMORIES value is not a number [%s]. Default value is 8. [bool - %s] " % (config[1], self.number(config[1])))
                    elif config[0] == "PROCSSHARED": # Procsshared
                        if self.number(config[1]) == True:
                            self.procsshared = int(config[1])
                        else:
                            log.write("Warning: PROCSSHARED value is not a number [%s]. Default value is 8. [bool - %s] " % (config[1], self.number(config[1])))
                    elif config[0] == "LEVELOFTHEORY": # leveloftheory for gaussian
                        self.leveloftheory = config[1] #no check due to the vast diversity
                    elif config[0] == "SPIN": # spin of the molecule
                        if self.number(config[1]) == True:
                            self.spin = int(config[1])
                        else:
                            log.write("Warning: SPIN value is not a number [%s]. Default value is 1. [bool - %s] " % (config[1], self.number(config[1])))
                    elif config[0] == "SOFTWARE": # SOFTWARE to use to optimise
                        self.software = config[1].upper()
                    elif config[0] == "RJCT": # RJCT is highly important. Define how close 2 atoms can be according to their diameter.
                        if self.number(config[1]) == True:
                            self.RJCT = float(config[1])
                        else:
                            log.write("Warning: RJCT value is not a number [%s]. Default value is 0.5 meaning the 2 atoms are supperimposed and the radii is going from one core to the other. [bool - %s] " % (config[1], self.number(config[1])))
                    elif config[0] == "TEMPERATURE": # Temperature for the Boltzmann distribution
                        if self.number(config[1]) == True:
                            self.Temperature = int(config[1])
                        else:
                            log.write("Warning: TEMPERATURE value is not a number [%s]. Default value is 298 K [bool - %s] " % (config[1], self.number(config[1])))
                    elif config[0] == "ENERGYWINDOW_CUTOFF": # Energy window cutoff for relevant conformers to be considered in the min, max, and median values for wB1, wB5 and wL
                        # ENERGYWINDOW_CUTOFF is a list. = xx xx xx.
                        config[1] = config[1].split()
                        for k in range(len(config[1])):
                            config[1][k] = config[1][k].strip(' ').strip('\t')
                            if self.number(config[1][k]) == True:
                                self.energywindow_cutoff.append(float(config[1][k]))
                            else: log.write("Warning: ENERGYWINDOW_CUTOFF value is a list of numbers in kcal/mol [submit: %s]." % config[1][k])
                    elif config[0] == "PRINT_CUTOFF": # Energy cutoff for relevant conformers to be included in log file. Useful when dealing with many conformers.
                        if self.number(config[1]) == True:
                            self.print_cutoff = float(config[1])
                        else:
                            log.write("Warning: PRINT_CUTOFF value is not a number [%s]. Default value is 5.0 Kcal/mol [bool - %s] " % (config[1], self.number(config[1])))
                    elif config[0] == "ANGLE_COUNT": # The division number to divide 360 degree in a certain number of time. angle_count = x; x calculated angles per dihedral. 360/x degrees = angle_step.
                        if self.number(config[1]) == True:
                            self.angle_count = int(config[1])
                        else:
                            log.write("Warning: ANGLE_COUNT value is not a number [%s]. Default value is 5 [bool - %s] " % (config[1], self.number(config[1])))
                    elif config[0] == "SPC_LEVELOFTHEORY": # Single Point Calculation level of theory after optimisation
                        self.singlepointcalculation = config[1] #no check due to the vast diversity
                    elif config[0] == "ATOMIC_MODEL": # atomic model for filter_gen and sterimol
                        radii = config[1].lower()
                        if radii != "cpk" and radii != "bondi":
                            radii = "cpk"
                            log.write("Warning: Radii must be \"cpk\" or \"bondi\". CPK atomic model is used by default.")
                        self.radii = config[1]
                    else:
                        log.write("Warning: Unknown setting in setup.ini. Check spelling [%s]" % config[0])
        else:
            self.loaded = False
            log.write("Warning: can't find setup.ini. Default values are used. [%s]" % self.path)

    def isLoaded(self):
        return self.loaded

    def number(self, s):
        try:
            float(s) # for int, long, float and complex
            return True
        except ValueError:
            return False
        return True
