#Pymol API
from pymol import cmd

########################################
########## w S T E R I M O L ###########
########################################

# Generate the wSterimol parameters from the optimised structures
# Use in Pymol command prompt:
# run wSterimol.py
# wSterimol dihedrals, atomid1, atomid2, (radii, walltime, directory, setup_path, verbose, force)

def wSterimol(dihedrals, atomid1 = 1, atomid2 = 2, directory = "temp", setup_path = '.', walltime = 300,  verbose = "False", force = False):
    # Log generation
    log_path = "log-%s" % (datetime.date.today())
    if os.path.exists(log_path+".pylog"):
        print("Warning: Continuing previous log file [%s.pylog]" % log_path)
    wlog = Log(log_path,"pylog")
    # Do weighted Sterimol
    generate(dihedrals, directory, setup_path, verbose, force)
    filter_gen(directory, setup_path, verbose)
    prepare_file(directory, setup_path, verbose)
    optimisation(directory, walltime, verbose, setup_path)
    filter_opt(directory, setup_path, verbose)
    Sterimol(atomid1, atomid2, directory, setup_path, verbose)
    weight(setup_path, verbose)
    wlog.write("---- wSterimol finished ----\n")
    wlog.finalize()
cmd.extend("wSterimol",wSterimol)
