#!/usr/bin/python
import subprocess, sys, os, math, datetime, shutil, getpass, datetime
from os import listdir
from os.path import isfile, join
from subprocess import Popen, PIPE

from pymol import cmd

########################################
########### F I L T E R - O ############
########################################

# Filter optimised structures generated from Mopac or Gaussian
# Use in Pymol command prompt:
# run log.py
# run ccParse.py
# run setup.py
# run filter_opt.py
# filter_opt (directory, verbose, setup_path)
# example: 

def filter_opt(directory = "temp", setup_path = "default", verbose = "False"):
    # If the directory exists
    if os.path.exists(directory):
        # Log generation
        log_path = "log-%s" % (datetime.date.today())
        if os.path.exists(log_path+".pylog"): 
            print "Warning: Continuing previous log file [%s.pylog]" % log_path
        log = Log(log_path,"pylog")
        log.write("\n\n########################################\n########### F I L T E R - O ############\n########################################\n\n")
        #verbose
        if verbose.lower() in ['true', '1', 't', 'y', 'yes']: verbose = True
        else: verbose = False
        # Retrieve all the files in the directory which are output with .log extension from gaussian or .out extension from mopac
        files = [f for f in listdir(directory) if isfile(join(directory, f)) and (f.split(".")[-1] == "log" or f.split(".")[-1] == "out") and len(f.split(".")) > 1 ]
        if len(files) > 0:
            # get setup.ini
            setup = Setup(log, setup_path)
            if setup.isLoaded() == True:
                # if optimisation, remove identical conformers
                if setup.scf != "1SCF":
                    log.write("\n\nRemove identical structure with an RMSD cutoff of: %s Angstroms " % setup.rmsd_cutoff)
                    #prepare bins
                    bins = [] #first of the bin is the median of the bin, ie the representation of the bin
                    energies = []
                    # Retrieve energy and create associated pdb
                    for i in range(len(files)):
                        fileData = getoutData(join(directory, files[i]))
                        #there is an energy, there was no error
                        if hasattr(fileData, "ENERGY"): 
                            export_pdb(directory, files[i], fileData)
                            energies.append(fileData.ENERGY) # hartree
                        else: #error: the optimisation probably failed for this conformer, we remove it and trigger a warning
                            os.remove(join(directory, files[i]))
                            log.write("WARNING: No energy was found. Removed from conformers [%s]" % join(directory, files[i]));
                    #reload
                    files = [f for f in listdir(directory) if isfile(join(directory, f)) and (f.split(".")[-1] == "log" or f.split(".")[-1] == "out") and len(f.split(".")) > 1 ]
                    # load all the pdb files
                    for file in files:
                        filesplit = file.split(".")
                        full_file_name = join(directory, "%s_OPT.pdb" % filesplit[0])
                        cmd.load( full_file_name, "%s_copy" % cmd.get_object_list()[0])
                    # calculate for each file
                    for i in range(len(files)):
                        log.write("\n%s_copy__%s RMSD:" % (cmd.get_object_list()[0],i+1), verbose)
                        RMSD = cmd.intra_fit("%s_copy"  % cmd.get_object_list()[0], i+1)
                        log.write("%s" % RMSD, verbose)
                        min_id = -1
                        min_rmsd = setup.rmsd_cutoff
                        for j in range(len(bins)):
                            current_rmsd = RMSD[bins[j][0]]
                            # if under a certain cut-off, considers the same structure
                            if current_rmsd < setup.rmsd_cutoff:
                                # put the structure in the bin with the minimum cutoff
                                if current_rmsd < min_rmsd:
                                    min_id = j
                                    min_rmsd = current_rmsd
                        # check results, and do accordingly
                        if min_id == -1: 
                            bins.append([i]) # doesn't belong to a bin, create a new bin
                            log.write("Create a new bin (#%s) for [%s]" % ((len(bins)-1), files[i]), verbose)
                        else: 
                            bins[min_id].append(i) #belongs to bin[min_id]
                            log.write("[%s] Added to bin #%s  with RMSD of %-6.4f Angstroms" % (files[i], min_id, min_rmsd), verbose)
                    #finally, filter by keeping one structure for each bin
                    for i in range(len(bins)):
                        min_energy = bins[i][0]
                        for j in range(1, len(bins[i])):
                            # keep the smallest energy conformation
                            if energies[bins[i][j]] < energies[min_energy]:
                                try: # remove previous files
                                    filename = files[min_energy].split(".")[0]
                                    # remove *.pdb, *.out or *.log
                                    os.remove(join(directory, "%s_OPT.pdb" % filename))
                                    os.remove(join(directory, files[min_energy]))
                                    if setup.software == "MOPAC": # remove *.mop, *.arc
                                        os.remove(join(directory, "%s.mop" % filename))
                                        os.remove(join(directory, "%s.arc" % filename))
                                    elif setup.software == "GAUSSIAN":
                                        os.remove(join(directory, "%s.com" % filename))
                                except: log.write("WARNING: Failed to remove [%s]" % join(directory, files[min_energy])); pass
                                min_energy = bins[i][j]
                            else:
                                try:  # remove these files
                                    filename = files[bins[i][j]].split(".")[0]
                                    # remove *.pdb, *.out or *.log
                                    os.remove(join(directory, "%s_OPT.pdb" % filename))
                                    os.remove(join(directory, files[bins[i][j]]))
                                    if setup.software == "MOPAC": # remove *.mop, *.arc
                                        os.remove(join(directory, "%s.mop" % filename))
                                        os.remove(join(directory, "%s.arc" % filename))
                                    elif setup.software == "GAUSSIAN":
                                        os.remove(join(directory, "%s.com" % filename))
                                except: log.write("WARNING: Failed to remove [%s]" % join(directory, files[bins[i][j]]));pass
                    #remove loaded structures from Pymol
                    cmd.delete("%s_copy" % cmd.get_object_list()[0])
                # create output
                output = open("energies.txt", 'w' )
                output.write("*******************************************************************************\n")
                output.write("** For non-commercial use only                                  Version %.2f **\n" % version)
                output.write("*******************************************************************************\n")
                output.write("** Cite this program as:  to be determined                                   **\n")
                output.write("** Paton Computational Chemistry Group,   web: http://paton.chem.ox.ac.uk    **\n")
                output.write("*******************************************************************************\n")
                output.write("**                                                                           **\n")
                output.write("**                             E N E R G I E S                               **\n")
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
                output.write("* Clusterisation cut-off was set at: %4.2f Angstroms                           *\n" % setup.rmsd_cutoff )
                output.write("* Angle count: %-3.f      | Atomic model: %-5s      | Temperature: %4.f K      *\n" % (setup.angle_count,setup.radii,setup.Temperature))
                output.write("*******************************************************************************\n")
                output.write(" File created on %s\n" % str(datetime.date.today()))
                output.write("\n Structure                       Etot (Hartree)")
                files = [f for f in listdir(directory) if isfile(join(directory, f)) and (f.split(".")[-1] == "log" or f.split(".")[-1] == "out") and len(f.split(".")) > 1 ]
                for filename in files:
                    out = getoutData(join(directory, filename ))
                    message = "\n %-30s" % filename +" %12.6f" % out.ENERGY
                    output.write(message)
                # terminate
                output.close()
                log.write("----------------------------\n---- Normal Termination ----\n----------------------------\n")
                log.finalize()
            else:
                log.write("Error: Failed to load setup.ini in [%s]. Fix it to continue." % setup_path)
        else: log.write("Error: No file to use in the directory [%s]" % directory)
    else: print "FATAL ERROR: Specified directory doesn't exist [%s]" % directory
        
cmd.extend("filter_opt",filter_opt)

# If .out files already exist, then we remove them to avoid conflict problems later.
def export_pdb(directory, file, Molspec):
    filesplit = file.split(".")
    full_file_name = join(directory, "%s_OPT.pdb" % filesplit[0])
    if os.path.exists(full_file_name): 
        os.remove(full_file_name)  # remove file
    file = open(full_file_name, 'w' )
    file.write("REMARK   1 File created by energy.py\n")
    for i in range(len(Molspec.CARTESIANS)):
        message = "HETATM %4.f  %s           0     % 7.3f % 7.3f % 7.3f                       %s\n" % (i+1, Molspec.ATOMTYPES[i], Molspec.CARTESIANS[i][0], Molspec.CARTESIANS[i][1], Molspec.CARTESIANS[i][2], Molspec.ATOMTYPES[i])
        file.write(message)
    file.write("END\n")
    file.close()
        