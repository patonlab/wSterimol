#!/usr/bin/python
from __future__ import print_function, absolute_import
import subprocess, sys, os, math, datetime, shutil, datetime
from os import listdir
from os.path import isfile, join
from subprocess import PIPE

from pymol import cmd
#run ../../../..\code\getRMSD.py
#run ../../../..\code\ccParse.py
#getRMSD conformers


def getRMSD(directory = "temp"):
    # If the directory exists
    if os.path.exists(directory):
        files = [f for f in listdir(directory) if isfile(join(directory, f)) and (f.split(".")[-1] == "log" or f.split(".")[-1] == "out") and len(f.split(".")) > 1 ]
        # load all the pdb files
        for file in files:
            filesplit = file.split(".")
            full_file_name = join(directory, "%s_OPT.pdb" % filesplit[0])
            cmd.load( full_file_name, "%s_copy" % cmd.get_object_list()[0])
        # calculate for the reference object
        RMSD = cmd.intra_fit("%s_copy"  % cmd.get_object_list()[0], 0)
        #print RMSD
        energies = []
        for i in range(len(files)):
            out = getoutData(join(directory, files[i] ))
            energies.append(out.ENERGY)
        # scale and convert hartree to J/mol. Then calculate Boltzman distribution
        hartree_to_Jmol = 2600 *1000 # hartree to J/mol
        Jmol_to_kcalmol = 0.24 /1000 # J/mol to kcal/mol
        min_energy = min(energies)
        for i in range(len(energies)):
            energies[i] = (energies[i] - min_energy) * hartree_to_Jmol*Jmol_to_kcalmol
        #write csv
        output = open("RMSD.csv", 'w' )
        output.write("Structure;Energy;RMSD\n")
        message = "%s;%.4f;%.4f\n" % (files[0], energies[0], 0) # RMSD = -1 for the same structure by default, but in our case 0 is more meningful
        output.write(message)
        for i in range(1,len(files)):
            message = "%s;%.4f;%.4f\n" % (files[i], energies[i], RMSD[i])
            output.write(message)
        # terminate
        output.close()
    else: print("FATAL ERROR: Specified directory doesn't exist [%s]" % directory)

cmd.extend("getRMSD",getRMSD)
