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
# generate [dihedral 1, .. ], (directory, setup_path, verbose)
# example: generate [[id 1, id 2, id 3, id 6],[id 2, id 3, id 6, id 11]], conformers

def generate(dihedrals, directory = "temp", setup_path = "default", verbose = "False"):
    # dihedrals = [[atom a, atom b, atom c, atom d]]
    # Log generation
    log = Log()
    log.write("\n\n########################################\n########### G E N E R A T E ############\n########################################\n\n")
    #verbose
    if verbose.lower() in ['true', '1', 't', 'y', 'yes']: verbose = True
    else: verbose = False
    # dihedrals: auto?
    if dihedrals.lower() in ['auto']: dihedrals = 'auto' #To DO: generate dihedrals automatically
    # get setup.ini
    setup = Setup(log, setup_path)
    if setup.isLoaded() == True:
        angle_step = 360 / setup.angle_count
        dihedrals = dihedrals.strip(' []').split(',')
        for i in range(0,len(dihedrals)):
            dihedrals[i]=dihedrals[i].strip(' []')
        if len(dihedrals)%4 != 0 and len(dihedrals) != 1:
            log.write("Error: Problem in dihedral selection. 4 selections are needed for each dihedral. [%s] " % dihedrals)
            return False
        dihedral_count=int(len(dihedrals)/4)
        combination_count=int(setup.angle_count**dihedral_count)
        log.write("Number of dihedrals: %s" % dihedral_count)
        log.write("For each dihedral: %s steps of %s degree" % (setup.angle_count, angle_step))
        log.write("That represents %s combinations" % combination_count)
        result = backup_generate(directory) # create backup
        if result: log.write("Making archives: %s" % result, verbose)
        else: log.write("Making directory: %s" % directory, verbose)
        if combination_count == 0:
            path = "%s/%s_0.pdb" % (directory, cmd.get_object_list()[0])
            try:
                cmd.save(path)
            except:
                log.write("Error: Can't save %s first conformer. One conformer will be missing!" % path, verbose)
        else:
            # Generate all the combinations
            combination_list = [{} for k in range(0,combination_count)]
            for i in range(0,combination_count): # iterate all the combinations of dihedrals
                dihedrals_curr_cbn=[{} for k in range(dihedral_count)]
                count = i
                for j in range(0,dihedral_count): # modify the dihedral according to the combination
                    dihedrals_curr_cbn[j]=int(count%setup.angle_count)
                    count=int(count/setup.angle_count)
                log.write("Combination %s : %s " % (i,dihedrals_curr_cbn), verbose)
                combination_list[i] = dihedrals_curr_cbn
            # Too many combinations?
            if combination_count > setup.combination_limit:
                todestroy_count = combination_count - setup.combination_limit
                log.write("The amount of possibility is %s. The limit is set at %s. %s combinations will be deleted to go down to the limit using %s algorithm." % (combination_count, setup.combination_limit, todestroy_count, setup.combination_remove_algorithm))
                if setup.combination_remove_algorithm == "random":
                    for i in range(0, todestroy_count): # loop as many time as the amount to randomly destroy
                        del combination_list[random.randint(0, len(combination_list))]
                log.write("%s combinations. Limit is %s. Reduction of combinations done." % (len(combination_list), setup.combination_limit), verbose)
            # Generate the files
            for i in range(len(combination_list)): # going through each combination
                for j in range(len(combination_list[i])): # going through each dihedral of the combination
                    try:
                        cmd.set_dihedral(dihedrals[j*4],dihedrals[j*4+1],dihedrals[j*4+2],dihedrals[j*4+3],combination_list[i][j]*angle_step)
                    except:
                        log.write("Error: Can't set the dihedrals. cmd.set_dihedrals failed.")
                log.write("Writing Combination %s : %s " % (i,combination_list[i]), verbose)
                path = "%s/%s_%s.pdb" % (directory, cmd.get_object_list()[0], str(i))
                try:
                    cmd.save(path)
                except:
                    log.write("Error: Can't save %s. One conformer will be missing!" % path)
                for j in range(0,dihedral_count):
                    cmd.undo()
        log.write("----------------------------\n---- Normal Termination ----\n----------------------------\n")
        log.finalize()
        return True
    else: 
        log.write("Error: Failed to load setup.ini in [%s]. Fix it to continue." % setup_path)
        return False

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
