#!/usr/bin/env python3

import subprocess
import os
import sys

###############
# FUNCTIONS ###
###############

def error(msg):
    sys.stderr.write("ERROR: {}\n".format(msg))
    sys.stderr.flush()
    sys.exit(1)

def execute(command):

    print('Executing {}'.format(command) + os.linesep)
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True)
    (stdout, stderr) = process.communicate()

    print(stderr.decode() + os.linesep, stdout.decode() + os.linesep)

execute('module load tacc-singularity') #&& singularity help
execute('singularity help')
