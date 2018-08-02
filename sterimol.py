#!/usr/bin/python

#Python Libraries 
import subprocess, sys, os
from os import listdir
from os.path import isfile, join
import datetime

#Pymol API
from pymol import cmd

########################################
########### S t E R I M O L ############
########################################

# Generate the Sterimol parameters from the optimised structures
# Use in Pymol command prompt:
# run log.py
# run ccParse.py
# run setup.py
# run sterimoltools.py
# run sterimol.py
# sterimol atomid1, atomid2, (directory, setup_path, verbose)

def Sterimol(atomid1 = 1, atomid2 = 2, directory = "temp", setup_path = "default", verbose = "False"):
    # If the directory exists
    if os.path.exists(directory):
        # Log generation
        log_path = "log-%s" % (datetime.date.today())
        if os.path.exists(log_path+".pylog"): 
            print "Warning: Continuing previous log file [%s.pylog]" % log_path
        log = Log(log_path,"pylog")
        log.write("\n\n########################################\n########### S t E R I M O L ############\n########################################\n\n")
        # Check arguments to avoid an error later
        if verbose.lower() in ['true', '1', 't', 'y', 'yes']:
            verbose = True
        else:
            verbose = False
        try:
            atomid1 = int(atomid1) # in Pymol, arguments are strings by default. Convert to int and check that is right.
            atomid2 = int(atomid2)
        except ValueError:
            log.write("FATAL ERROR: Atom id must be an integer. In Pymol, Label (L) -> atom identifiers -> ID\n")
            return
        # Retrieve all the files in the directory
        files = [f for f in listdir(directory) if isfile(join(directory, f)) and (f.split(".")[-1] == "log" or f.split(".")[-1] == "out") and len(f.split(".")) > 1 ]
        if len(files) > 0:
            # get setup.ini
            setup = Setup(log, setup_path)
            if setup.isLoaded() == True:
                output = open("Sterimol.txt", 'w' )
                output.write("*******************************************************************************\n")
                output.write("** For non-commercial use only                                  Version %.2f **\n" % version)
                output.write("*******************************************************************************\n")
                output.write("** Cite this program as:  to be determined                                   **\n")
                output.write("** Paton Computational Chemistry Group,   web: http://paton.chem.ox.ac.uk    **\n")
                output.write("*******************************************************************************\n")
                output.write("**                                                                           **\n")
                output.write("**                             S T E R I M O L                               **\n")
                output.write("**                                                                           **\n")
                output.write("*******************************************************************************\n")
                output.write("* sterimol -a1 %-3s -a2 %-3s -radii %-6s                                      *\n" % (atomid1, atomid2, setup.radii) )
                output.write("* Clusterisation cut-off was set at: %4.2f Angstroms                           *\n" % setup.rmsd_cutoff )
                output.write("* Angle count: %-3.f      | Atomic model: %-5s      | Temperature: %4.f K      *\n" % (setup.angle_count,setup.radii,setup.Temperature))
                output.write("*******************************************************************************\n")
                output.write(" File created on %s | A for Angstroms\n\n" % str(datetime.date.today()))
                message = " Structure                          L1 (A)   B1 (A)   B5 (A) \n"
                output.write(message)
                log.write(message, verbose)
                for filename in files:
                    filesplit = filename.split(".")
                    #prepare the job
                    try:
                        file_Params = calcSterimol(join(directory, filename), setup.radii, atomid1, atomid2, verbose)
                    except ValueError:
                        log.write("FATAL ERROR: An error occured in Sterimol calculation.\n\n%s\n\n" % ValueError)
                        return
                    lval = file_Params.lval; B1 = file_Params.B1; B5 = file_Params.newB5
                    message = " %-31s " % filename+" %8.2f" % lval+ " %8.2f" % B1+ " %8.2f" % B5
                    log.write(message, verbose)
                    output.write("%s\n" % message)
                output.close()
                log.write("----------------------------\n---- Normal Termination ----\n----------------------------\n")
                log.finalize()
            else: log.write("Error: Failed to load setup.ini in [%s]. Fix it to continue." % setup_path)
        else: log.write("Error: No file in the directory [%s]" % directory)
    else: print "FATAL ERROR: Specified directory doesn't exist [%s]" % directory
        
cmd.extend("sterimol",Sterimol)


