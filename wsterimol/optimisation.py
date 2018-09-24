#!/usr/bin/python
from __future__ import print_function, absolute_import
import subprocess, sys, os, math, datetime, shutil, platform
from os import listdir
from os.path import isfile, join
from subprocess import PIPE


from pymol import cmd

########################################
######## O P T M I S A T I O N #########
########################################

# Optimisation of the structures via Pymol with MOPAC
# Use in Pymol command prompt:
# run wSterimol.py
# run sterimoltools.py
# run setup.py
# run optimisation.py
# optimisation (directory, walltime, verbose, setup_path)
# example: optimisation conformers

def optimisation(directory = "temp", walltime = 300, verbose = "False", setup_path = "default", silentmode = "True"):
    # If the directory exists
    if os.path.exists(directory):
        # Log generation
        log_path = "log-%s" % (datetime.date.today())
        if os.path.exists(log_path+".pylog"):
            print("Warning: Continuing previous log file [%s.pylog]" % log_path)
        log = Log(log_path,"pylog")
        log.write("\n\n########################################\n###########   E N E R G Y   ############\n########################################\n\n")
        #verbose
        if verbose.lower() in ['true', '1', 't', 'y', 'yes']: verbose = True
        else: verbose = False
        #silentmode
        if silentmode.lower() in ['true', '1', 't', 'y', 'yes']: silentmode = True
        else: silentmode = False
        # Retrieve all the files in the directory
        files = [f for f in listdir(directory) if isfile(join(directory, f))]
        if len(files) > 0:
            # get setup.ini
            setup = Setup(log, setup_path)
            if setup.isLoaded() == True:
                # check walltime is ok
                try:
                    cycles = int(int(walltime)/0.3) # sleep time is 0.3 sec, so 1000 sec is 1000 sec / 0.3 sec = number of cycles
                except ValueError:
                    log.write("Error: Walltime keyword is not a number. Default value is set at 300 seconds [%s]" % walltime)
                    cycles = 1000 # walltime: 300s by default, 1000 cycles
                if setup.software == "MOPAC": # MOPAC software
                    # calculate for each file
                    for filename in files:
                        filesplit = filename.split(".")
                        if len(filesplit) > 1: # there is an extension
                            if filesplit[-1] == "mop":
                                # check output file doesn't already exist
                                if backup_energy(directory, filename, log) == True:
                                    log.write("\n-------------------------\nJob starts [%s]" % join(directory, filename), verbose)
                                    # calculate the energy
                                    startupinfo = None
                                    if platform.system() == "Windows" and silentmode == True:
                                        # Prevent .exe window from popping up in Windows
                                        startupinfo = subprocess.STARTUPINFO()
                                        startupinfo.dwFlags |= subprocess.STARTF_USESTDHANDLES | subprocess.STARTF_USESHOWWINDOW
                                        startupinfo.wShowWindow = subprocess.SW_HIDE
                                    job = subprocess.call([setup.exe,  join(directory, filename)], stdout=PIPE, startupinfo=startupinfo)
                                    # Make sure optimization is complete
                                    n = 0
                                    while(isJobFinished(directory, filename) == False and n < cycles):
                                        time.sleep(0.3)
                                        n = n +1
                                    if isJobFinished(directory, filename) == True:
                                        filesplit = filename.split(".")
                                        out = getoutData(join(directory, "%s.out" % filesplit[0]))
                                        message = "   %-30s" % (filesplit[0]+".out") +"  %12.6f" % out.ENERGY
                                        log.write("Job Finished [%s]" % join(directory, filename), verbose)
                                        log.write("   Structure                        Etot (Hartree)", verbose)
                                        log.write(message, verbose)
                                        log.write("-------------------------", verbose)
                                    else:
                                        log.write("Warning: Exit abnormally. Look at the output file of the input [%s]" % join(directory, filename))
                                else:
                                    log.write("Warning: Failed to delete the existing output of [%s]. Skip Job calculation!" % join(directory, filename))
                    # terminate
                    log.write("----------------------------\n---- Normal Termination ----\n----------------------------\n")
                    log.finalize()
                elif setup.software == "GAUSSIAN": # MOPAC software
                    # calculate for each file
                    for filename in files:
                        filesplit = filename.split(".")
                        if len(filesplit) > 1: # there is an extension
                            if filesplit[-1] == "com":
                                # check output file doesn't already exist
                                #if backup_energy(directory, filename, log) == True:
                                    log.write("\n-------------------------\nJob starts [%s]" % join(directory, filename), verbose)
                                    # calculate the energy
                                    startupinfo = None
                                    if platform.system() == "Windows" and silentmode == True:
                                        # Prevent .exe window from popping up in Windows
                                        startupinfo = subprocess.STARTUPINFO()
                                        startupinfo.dwFlags |= subprocess.STARTF_USESTDHANDLES | subprocess.STARTF_USESHOWWINDOW
                                        startupinfo.wShowWindow = subprocess.SW_HIDE

                                    command = setup.exe + ' ' + join(directory, filename.split(".")[0])
                                    #log.write(command)
                                    job = subprocess.call(command, stdout=PIPE, shell=True)
                                    # Make sure optimization is complete
                                    n = 0
                                    #while(isJobFinished(directory, filename) == False and n < cycles):
                                    #    time.sleep(0.3)
                                    #    n = n +1
                                    if isJobFinished(directory, filename) == True:
                                        filesplit = filename.split(".")
                                        out = getoutData(join(directory, "%s.out" % filesplit[0]))
                                        message = "   %-30s" % (filesplit[0]+".out") +"  %12.6f" % out.ENERGY
                                        log.write("Job Finished [%s]" % join(directory, filename), verbose)
                                        log.write("   Structure                        Etot (Hartree)", verbose)
                                        log.write(message, verbose)
                                        log.write("-------------------------", verbose)
                                    else:
                                        log.write("Warning: Exit abnormally. Look at the output file of the input [%s]" % join(directory, filename))
                                #else:
                                #    log.write("Warning: Failed to delete the existing output of [%s]. Skip Job calculation!" % join(directory, filename))
                    # terminate
                    log.write("----------------------------\n---- Normal Termination ----\n----------------------------\n")
                    log.finalize()
                else:
                    log.write("Error: Wanted software is not compatible with this energy calculation script. Use another software independantly. [%s]" % setup.software)
            else:
                log.write("Error: Failed to load setup.ini in [%s]. Fix it to continue." % setup_path)
        else: log.write("Error: No file to use in the directory [%s]" % directory)
    else: print("FATAL ERROR: Specified directory doesn't exist [%s]" % directory)


cmd.extend("optimisation",optimisation)

# Check that a computational chemistry job has finished #######
def isJobFinished(directory, file):
    filesplit = file.split(".")
    if not os.path.exists(join(directory, "%s.out" % filesplit[0])):
        return False
    else:
        outfile = open(join(directory, "%s.out" % filesplit[0]),"r")
        inlines = outfile.readlines()
        outfile.close()
        for line in inlines:
            if line.find("== MOPAC DONE ==") > -1: return True
            if line.find("Normal termination of Gaussian") > -1: return True
        return False

# If .out files already exist, then we remove them to avoid conflict problems later.
def backup_energy(directory, file, log):
    filesplit = file.split(".")
    if not os.path.exists(join(directory, "%s.out" % filesplit[0])):
        return True #file doesn't exist. do nothing
    else: # file exist, we must delete it
        new_path = "%s/%s" % (directory, datetime.date.today())  # archive directory
        if not os.path.exists(new_path):  # create archive folder if it does not exist
            os.makedirs(new_path)

        full_file_name = join(directory, "%s.out" % filesplit[0])
        try:
            shutil.copy(full_file_name, new_path)  # archive file
        except:
            log.write("Error: Failed to copy file [%s] to [%s]" % (full_file_name, new_path))
            pass
        try:
            os.remove(full_file_name)  # remove file if copied successfully or not
        except:
            return False
        return True
