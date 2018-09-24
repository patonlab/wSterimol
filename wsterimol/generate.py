#!/usr/bin/python
from __future__ import print_function, absolute_import

import subprocess, sys, os, math, datetime, shutil

from pymol import cmd

########################################
########### G E N E R A T E ############
########################################

# Generate possible conformer structures via Pymol
# Use in Pymol command prompt:
# run wSterimol.py
# run setup.py
# run generate.py
# generate [dihedral 1, .. ], (directory, setup_path, verbose, force)
# example: generate [[id 1, id 2, id 3, id 6],[id 2, id 3, id 6, id 11]], conformers

def generate(dihedrals, directory = "temp", setup_path = "default", verbose = "False", force = False):
    # dihedrals = [[atom a, atom b, atom c, atom d]]
    # Log generation
    log_path = "log-%s" % (datetime.date.today())
    if os.path.exists(log_path+".pylog"):
        print("Warning: Continuing previous log file [%s.pylog]" % log_path)
    log = Log(log_path,"pylog")
    log.write("\n\n########################################\n########### G E N E R A T E ############\n########################################\n\n")
    #verbose
    if verbose.lower() in ['true', '1', 't', 'y', 'yes']: verbose = True
    else: verbose = False
    # get setup.ini
    setup = Setup(log, setup_path)
    if setup.isLoaded() == True:
        angle_step = 360 / setup.angle_count
        dihedrals = dihedrals.strip(' []').split(',')
        for i in range(0,len(dihedrals)):
            dihedrals[i]=dihedrals[i].strip(' []')
        if len(dihedrals)%4 != 0 and len(dihedrals) != 1:
            log.write("Error: Problem in dihedral selection. 4 selections are needed for each dihedral. [%s] " % dihedrals)
            return
        dihedral_count=int(len(dihedrals)/4)
        conbinations=int(setup.angle_count**dihedral_count)
        log.write("Number of dihedrals: %s" % dihedral_count)
        log.write("For each dihedral: %s steps of %s degree" % (setup.angle_count, angle_step))
        log.write("That represents %s combinations" % conbinations)
        if conbinations >= 1000:
            log.write("You are about to generate %s files. The calculation has been canceled.\nForce the system if you want to continue:\ngenerate [dihedral 1, .. ], (directory, setup_path, verbose, force)" % conbinations)
            if not force:
                return
            else: log.write("You are about to generate %s files. System has been forced to continue the calculation." % conbinations)
        result = backup_generate(directory) # create backup
        if result: log.write("Making archives: %s" % result, verbose)
        else: log.write("Making directory: %s" % directory, verbose)
        if conbinations == 0:
            path = "%s/%s_0.pdb" % (directory, cmd.get_object_list()[0])
            try:
                cmd.save(path)
            except:
                log.write("Error: Can't save %s first conformer. One conformer will be missing!" % path, verbose)
        else:
            for i in range(0,conbinations): # iterate all the combinations of dihedrals
                dihedrals_curr_cbn=[{} for k in range(dihedral_count)]
                if i >= 0:
                    count = i
                    for j in range(0,dihedral_count): # modify the dihedral according to the conbination
                        dihedrals_curr_cbn[j]=int(count%setup.angle_count)
                        count=int(count/setup.angle_count)
                        try:
                            cmd.set_dihedral(dihedrals[j*4],dihedrals[j*4+1],dihedrals[j*4+2],dihedrals[j*4+3],dihedrals_curr_cbn[j]*angle_step)
                        except:
                            log.write("Error: Can't set the dihedrals. cmd.set_dihedrals failed.")
                log.write("Combination %s : %s " % (i,dihedrals_curr_cbn), verbose)
                path = "%s/%s_%s.pdb" % (directory, cmd.get_object_list()[0], str(i))
                try:
                    cmd.save(path)
                except:
                    log.write("Error: Can't save %s. One conformer will be missing!" % path)
                if i >= 0: # undo the dihedrals change for each dihedral before doing the next combination
                    for j in range(0,dihedral_count):
                        cmd.undo()
        log.write("----------------------------\n---- Normal Termination ----\n----------------------------\n")
        log.finalize()
    else: log.write("Error: Failed to load setup.ini in [%s]. Fix it to continue." % setup_path)

def backup_generate(path):
    new_path = "%s/%s" % (path, datetime.date.today())  # archive directory

    if os.path.exists(path):  # archive files in backup directory
        dir_files = [backup for backup in os.listdir(path) if os.path.isfile(os.path.join(path, backup))]  # existing files
        if len(dir_files) > 0:
            if not os.path.exists(new_path):  # create archive folder if it does not exist
                os.makedirs(new_path)

            for file_name in dir_files:
                full_file_name = os.path.join(path, file_name)
                try:
                    shutil.copy(full_file_name, new_path)  # archive file
                    os.remove(full_file_name)  # remove file if copied successfully
                except:
                    pass
            return new_path # backed up
    else:
        os.makedirs(path)  # create new directory if it doesn't exist
        return # created the directory

cmd.extend("generate",generate)
