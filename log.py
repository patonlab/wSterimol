#!/usr/bin/python
from __future__ import print_function, absolute_import

import sys, os
import timeit
from datetime import timedelta

# Weighted Sterimol Program version
version = 1.03

class Log:
    def __init__(self,filein,suffix):
        # Create the log file at the input path
        self.start_time = timeit.default_timer()
        self.log = open(filein+"."+suffix, 'a+' )

    # Write a message to the log and the terminal by default
    def write(self, message, verbose = True):
        # Print the message
        if verbose: print(message) # versatile, depend on the called function
        # Write to log
        self.log.write(message + "\n")

    # Write a fatal error, finalize and terminate the program
    def fatal(self, message):
        # Print the message
        print(message+"\n")
        # Write to log
        self.log.write(message + "\n")
        # Finalize the log
        self.finalize()
        # End the program
        sys.exit(1)

    # Finalize the log file
    def finalize(self):
        self.end_time = timeit.default_timer()
        self.log.write("Time elapsed during script execution: %s \n" % timedelta(seconds=self.end_time - self.start_time))
        self.log.close()
