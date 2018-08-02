#!/usr/bin/python
import subprocess, sys, os, math, datetime, shutil, getpass
from os.path import isfile, join

from pymol import cmd

########################################
###########   W E I G H T   ############
########################################

# Retrieve the energies and the Sterimol parameters in order to get the Boltzmann distribution
# Use in Pymol command prompt:
# run weight.py
# weight (verbose, setup_path)

def weight(setup_path = "default", verbose = "False"):
    # If files exist
    if os.path.exists("energies.txt") and os.path.exists("Sterimol.txt"):
        log_path = "log-%s" % (datetime.date.today())
        if os.path.exists(log_path+".pylog"): 
            print "Warning: Continuing previous log file [%s.pylog]" % log_path
        log = Log(log_path,"pylog")
        log.write("\n\n########################################\n###########   W E I G H T   ############\n########################################\n\n")
        #verbose
        if verbose.lower() in ['true', '1', 't', 'y', 'yes']: verbose = True
        else: verbose = False
        # get setup.ini
        setup = Setup(log, setup_path)
        if setup.isLoaded() == True:
            # create lists
            filenames_val = []
            energies_val = []
            sterimol_val = []
            # Retrieve energies 
            outenergy = open("energies.txt","r")
            outlinesenergy = outenergy.readlines()
            outenergy.close()
            start = -1
            for i in range(0,len(outlinesenergy)):
                if outlinesenergy[i].find("Structure") > -1: 
                    start = i+1
                    break
            if start >= 0:
                for i in range(start,len(outlinesenergy)):
                    split = outlinesenergy[i].split()
                    if len(split) >= 2: 
                        if is_str(split[0]):
                            if try_number(split[1]):
                                filenames_val.append(str(split[0]))
                                energies_val.append(float(split[1]))
                                sterimol_val.append([]) # fill to have the same number of rows
                            else: log.write("Error: energy value is not a number [%s]" % split[1])
                        else: log.write("Error: energy file name is not a string. Corrupted file.")
                    else: log.write("Error: energy values are corrupted. Not enough data.")
                if len(filenames_val) > 0 and len(energies_val) > 0 :
                    # Retrieve Sterimol
                    outsterimol = open("Sterimol.txt","r")
                    outlinessterimol = outsterimol.readlines()
                    outsterimol.close()
                    start = -1
                    for i in range(0,len(outlinessterimol)):
                        if outlinessterimol[i].find("Structure") > -1: 
                            start = i+1
                            break
                    if start >= 0:
                        for i in range(start,len(outlinessterimol)):
                            split = outlinessterimol[i].split()
                            if len(split) >= 4: 
                                if is_str(split[0]):
                                    if try_number(split[1]) and try_number(split[2]) and try_number(split[3]):
                                        for j in range(len(filenames_val)):
                                            if filenames_val[j] == str(split[0]): # find the same file name
                                                sterimol_val[j] = [ float(split[1]), float(split[2]), float(split[3])]
                                    else: log.write("Error: Sterimol values are not numbers")
                                else: log.write("Error: Sterimol file name is not a string. Corrupted file.")
                            else: log.write("Error: Sterimol values are corrupted. Not enough data.")
                    else: log.write("Error: cant find the start of Sterimol values")
                else: log.write("Error: No data was retrieved from energy values")
            else: log.write("Error: cant find the start of energy values")
            # calculate weighted sterimol
            if len(filenames_val) > 0 and len(energies_val) > 0 and len(sterimol_val) > 0:
                # check sterimol_val has been correctly filled up
                loaded = True
                for i in range(len(sterimol_val)):
                    if not len(sterimol_val[i]) == 3:
                        loaded = False
                        log.write("Error: Sterimol values are corrupted or missing [%s values] for the corresponding file [%s]." % (len(sterimol_val[i]), filenames_val[i]))
                if loaded == True:
                    distrib_val = []
                    R = 8.314 # J K−1 mol−1
                    hartree_to_Jmol = 2600 *1000 # hartree to J/mol
                    Jmol_to_kcalmol = 0.24 /1000 # J/mol to kcal/mol
                    min_energy = min(energies_val)
                    # scale and convert hartree to J/mol. Then calculate Boltzman distribution
                    for i in range(len(energies_val)):
                        energies_val[i] = (energies_val[i] - min_energy) * hartree_to_Jmol
                        distrib_val.append(round(math.exp(-energies_val[i]/(R* setup.Temperature ))*10000)/10000) # remove non significant numbers
                    # scale distribution. Somme(distributions) = 100%
                    sum_distrib = sum(distrib_val)
                    for i in range(len(distrib_val)):
                        distrib_val[i] = distrib_val[i] / sum_distrib
                    # Print data
                    output = open("weighted.txt", 'w' )
                    output.write("*******************************************************************************\n")
                    output.write("** For non-commercial use only                                  Version %.2f **\n" % version)
                    output.write("*******************************************************************************\n")
                    output.write("** Cite this program as:  to be determined                                   **\n")
                    output.write("** Paton Computational Chemistry Group,   web: http://paton.chem.ox.ac.uk    **\n")
                    output.write("*******************************************************************************\n")
                    output.write("**                                                                           **\n")
                    output.write("**                    W E I G H T E D   S T E R I M O L                      **\n")
                    output.write("**                                                                           **\n")
                    output.write("*******************************************************************************\n")
                    output.write("* CALCULATION WITH %-30s                             *\n" % setup.software)
                    if setup.software == "MOPAC": 
                        output.write("* MOPAC EXECUTIVE PATH: %-53s *\n" % setup.exe)
                        output.write("* using force-field: %-56s *\n" % setup.SE)
                        keywords = "%s charge=%s %s" % (setup.SE, setup.charge, setup.scf)
                        output.write("* %-75s *\n" % keywords)
                    elif setup.software == "GAUSSIAN": 
                        output.write("* %%mem=%2.fGB                                                                   *\n" % setup.memories )
                        output.write("* %%nprocshared=%2.f                                                             *\n" % setup.procsshared)
                        output.write("* # opt=(maxcycles=160) freq=noraman %-40s *\n" % setup.leveloftheory) 
                        output.write("* %2.f %2.f                                                                       *\n" % (setup.charge, setup.spin))
                        if setup.singlepointcalculation != "":
                            output.write("* --Link1--                                                                   *\n")
                            output.write("* # geom=check guess=read %-51s *\n" % setup.singlepointcalculation)
                            output.write("* %2.f %2.f                                                                       *\n" % (setup.charge, setup.spin))
                    output.write("*******************************************************************************\n")
                    output.write("* Clusterisation cut-off was set at: %4.2f Angstroms                           *\n" % setup.rmsd_cutoff )
                    output.write("* Angle count: %-3.f      | Atomic model: %-5s      | Temperature: %4.f K      *\n" % (setup.angle_count,setup.radii,setup.Temperature))
                    output.write("* Print cut-off: %5.1f kcal/mol                                               *\n" % (setup.print_cutoff))
                    for i in range(len(setup.energywindow_cutoff)):
                        output.write("* Energy window cut-off #%-2.f: %3.1f kcal/mol                                     *\n" % (i, setup.energywindow_cutoff[i]))
                    output.write("*******************************************************************************\n")
                    output.write(" File created on %s | A for Angstroms\n\n" % str(datetime.date.today()))
                    output.write(" Structures                   E (kcal/mol)  L1 (A)  B1 (A)  B5 (A)    (%)")
                    for i in range(len(filenames_val)):
                        if (energies_val[i] * Jmol_to_kcalmol) < setup.print_cutoff:
                            output.write("\n %-25s %15.2f %7.2f %7.2f %7.2f %6.2f" % (filenames_val[i], (energies_val[i] * Jmol_to_kcalmol), sterimol_val[i][0], sterimol_val[i][1], sterimol_val[i][2], (distrib_val[i]*100) ) )
                    # calculate final values
                    L1 = []; B1 = []; B5 = []
                    for i in range(len(sterimol_val)): # rearrange the lists differently to easily sum everything
                        L1.append(sterimol_val[i][0])
                        B1.append(sterimol_val[i][1])
                        B5.append(sterimol_val[i][2])
                    L1_weighted = 0
                    for i in range(len(L1)):
                        L1_weighted = L1_weighted + L1[i] * distrib_val[i]
                    B1_weighted = 0
                    for i in range(len(B1)):
                        B1_weighted = B1_weighted + B1[i] * distrib_val[i]
                    B5_weighted = 0
                    for i in range(len(B5)):
                        B5_weighted = B5_weighted + B5[i] * distrib_val[i]
                    output.write("\n\n*******************************************************************************\n")
                    output.write("**                                                                           **\n")
                    output.write("**                        wL (A)   wB1 (A)   wB5 (A)                         **\n")
                    output.write("**                       %7.2f   %7.2f   %7.2f                         **\n" % (L1_weighted,B1_weighted,B5_weighted))
                    output.write("**                                                                           **\n")
                    output.write("*******************************************************************************\n")
                    # calculate min max for each different energy windows
                    for k in range(len(setup.energywindow_cutoff)):
                        output.write("\n\n*******************************************************************************\n")
                        output.write("* Energy window cut-off #%-2.f: %3.1f kcal/mol                                     *\n" % (k, setup.energywindow_cutoff[k]))
                        #remove the very high values
                        filenames_val_copy = []
                        energies_val_copy = []
                        sterimol_val_copy = []
                        distrib_val_copy = []
                        for i in range(len(filenames_val)):
                            if (energies_val[i] * Jmol_to_kcalmol) < setup.energywindow_cutoff[k]:
                                filenames_val_copy.append(filenames_val[i])
                                energies_val_copy.append(energies_val[i])
                                sterimol_val_copy.append(sterimol_val[i])
                                distrib_val_copy.append(distrib_val[i])
                        # calculate final values
                        L1_copy = []
                        B1_copy = []
                        B5_copy = []
                        for i in range(len(sterimol_val_copy)): # rearrange the lists differently to easily sum everything
                            L1_copy.append(sterimol_val_copy[i][0])
                            B1_copy.append(sterimol_val_copy[i][1])
                            B5_copy.append(sterimol_val_copy[i][2])
                        output.write("*                                                                             *\n")
                        output.write("*        Lmin (A)  Lmax (A)  B1min (A)  B1max (A)  B5min (A)  B5max (A)       *\n")
                        output.write("*        %8.2f  %8.2f  %9.2f  %9.2f  %9.2f  %9.2f       *\n" % (min(L1_copy), max(L1_copy), min(B1_copy), max(B1_copy), min(B5_copy), max(B5_copy)))
                        output.write("*                                                                             *\n")
                        output.write("*******************************************************************************\n")
                        #show values
                        output.write(" List of the selected conformers accounting for %.3f%% of the conformers:\n" % (sum(distrib_val_copy)*100))
                        output.write(" Structures                   E (kcal/mol)  L1 (A)  B1 (A)  B5 (A)    (%)")
                        for i in range(len(filenames_val_copy)):
                            output.write("\n %-25s %15.2f %7.2f %7.2f %7.2f %6.2f" % (filenames_val_copy[i], (energies_val_copy[i] * Jmol_to_kcalmol), sterimol_val_copy[i][0], sterimol_val_copy[i][1], sterimol_val_copy[i][2], (distrib_val_copy[i]*100) ) )
                    # merge energies.txt and sterimol.txt then delete them
                    output.write("\n\n")
                    output.write("                         COPY OF DATA FOR ARCHIVES\n\n")
                    for i in range(0,len(outlinesenergy)):
                        output.write(outlinesenergy[i])
                    try: os.remove("energies.txt")
                    except NameError: log.write("WARNING: Failed to remove energies.txt.\n %s" % NameError)
                    output.write("\n\n")
                    for i in range(0,len(outlinessterimol)):
                        output.write(outlinessterimol[i])
                    try: os.remove("Sterimol.txt")
                    except NameError: log.write("WARNING: Failed to remove Sterimol.txt. %s" % NameError)
                    # exit
                    output.close()
                    log.write("\n----------------------------\n---- Normal Termination ----\n----------------------------\n")
                    log.finalize()
            else: log.write("Error: No data was retrieved for the weighted sterimol. Exit.")
        else: log.write("Error: Failed to load setup.ini in [%s]. Fix it to continue." % setup_path)
    else: print "FATAL ERROR: Files containing energies and Sterimol values don't exist"
        
   
cmd.extend("weight",weight)


def try_number(s):
    try:
        complex(s) # for int, long, float and complex
    except ValueError:
        return False
    return True

def is_str(s):
    try:
        str(s) # for str
    except ValueError:
        return False
    return True
