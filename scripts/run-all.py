#!/usr/bin/env python3

###############################################################################
#                                                                             #
#  run-all.py                                                                 #
#                                                                             #
#  A master script to run all the steps of radcot in sequence on slurm/uahpc  #
#                                                                             #
#    Copyright (C) Scott G Daniel                                             #
#                                                                             #
###############################################################################
#                                                                             #
#    This library is free software; you can redistribute it and/or            #
#    modify it under the terms of the GNU Lesser General Public               #
#    License as published by the Free Software Foundation; either             #
#    version 3.0 of the License, or (at your option) any later version.       #
#                                                                             #
#    This library is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU        #
#    Lesser General Public License for more details.                          #
#                                                                             #
#    You should have received a copy of the GNU Lesser General Public         #
#    License along with this library.                                         #
#                                                                             #
###############################################################################

__author__ = "Scott G Daniel"
__copyright__ = "Copyright 2018"
__credits__ = "Scott G Daniel"
__license__ = "LGPLv3"
__maintainer__ = "Scott G Daniel"
__email__ = "scottdaniel@email.arizona.edu"
__status__ = "Development"

"""Boilerplate stuff
This script essentially just parses out the metadata file
and the options files
and feeds them to each of the steps
"""

import os, sys, argparse
import plumbum
import pandas as pd

#WORK env var will be present on TACC
#But may not be set when testing locally
if os.getenv('WORK') is None:
    os.environ['WORK'] = './'

####################
# ARGUMENTS ########
####################
parser = argparse.ArgumentParser(description=
        "A master script to run all the other scripts\n"
        "To see what options are available, go to the scripts\n"
        "directory within each step and type [script name] -h\n"
        "Ex: under 01-centrifuge-patric/script, you can type:\n"
        "cfuge_to_patric.py to get the options for downloading\n"
        "PATRIC bacterial genomes.\n",
        formatter_class=argparse.RawTextHelpFormatter)

inputs = parser.add_argument_group('Required Inputs and Parameters')

inputs.add_argument('-m', '--metadata', 
        dest='metadata', metavar='FILENAME',
        default='metadata.txt',
        help="File containing file / sample information.\n"
        "Use the metadata_template to start and DO NOT change\n"
        "headings. [ Default = metadata.txt ]")

step_one = parser.add_argument_group('Step one option\'s files \(01-centrifuge-patric\)')

step_one.add_argument('-C', '--centrifuge-opts', 
        dest='cent_opts', metavar='FILENAME', 
        default=os.path.join(os.getenv('WORK'),'centrifuge-opts.txt'),
        help="File with additional options for centrifuge\n"
        "Please format your options like so:\n"
        "-o1 option1\n"
        "-o2 option2\n"
        "[ Default = $WORK/centrifuge-opts.txt ]")

step_one.add_argument('-P', '--PATRIC-opts', 
        dest='patric_opts', metavar='FILENAME', 
        default=os.path.join(os.getenv('WORK'),'PATRIC-opts.txt'),
        help="File with additional options for getting genomes from patric\n"
        "Please format your options like so:\n"
        "-o1 option1\n"
        "-o2 option2\n"
        "[ Default = $WORK/PATRIC-opts.txt ]")

step_two = parser.add_argument_group(
        'Step two option\'s files \(02-bowtie-samtools\)')

step_two.add_argument('-B', '--bowtie2-opts', 
        dest='bowtie2_opts', metavar='FILENAME',
        default=os.path.join(os.getenv('WORK'),'bowtie2-opts.txt'),
        help="File with additional options for bowtie2 alignment for RNA\n"
        "Please format your options like so:\n"
        "-o1 option1\n"
        "-o2 option2\n"
        "[ Default = $WORK/bowtie2-opts.txt ]")

step_three = parser.add_argument_group(
        'Step three option\'s files \(03-count-deseq\)')

step_three.add_argument('-H', '--htseq-count-opts', 
        dest='htseq_count_opts', metavar='FILENAME',
        default=os.path.join(os.getenv('WORK'),'htseq-opts.txt'),
        help="File with additional opts for htseq-count\n"
        "Options must be one per line like so:\n"
        "-o1 option1\n"
        "-o2 option2\n"
        "[ Default = $WORK/htseq-opts.txt ]")

step_three.add_argument('-D', '--deseq2-opts', 
        dest='deseq2_opts', metavar='FILENAME',
        default=os.path.join(os.getenv('WORK'),'deseq2-opts.txt'),
        help="File with additional opts for deseq2\n"
        "Options must be one per line like so:\n"
        "-o1 option1\n"
        "-o2 option2\n"
        "[ Default = $WORK/deseq2-opts.txt ]")

gen_opts = parser.add_argument_group('General Options') 

gen_opts.add_argument('-d', '--debug', action='store_true',
        help="Extra logging messages.") 

gen_opts.add_argument('-t', '--threads', 
        dest='threads', metavar='INT', 
        type=int, default=1,
        help="number of alignment threads to launch.\n"
        " [ Default = 1 ] ")

args = parser.parse_args()

#######################
# GENERAL FUNCTIONS ###
#######################

#really basic checker, check that options file exists and then check that each line begins with a '-', then parse
def parse_options_text(options_txt_path):
    if not (os.path.isfile(options_txt_path)):
        print("Options text {} does not exist or is not a file\n".format(options_txt_path))
        return None
    else:
        options_string = ''
        with open(options_txt_path) as options_txt:
            for line in options_txt:
                if not line.startswith('-'):
                    print("Skipping line that doesnt have hyphen\n")
                else:
                    options_string += line.replace('\n', ' ')

    if args.debug:
        print("These are the options for {}:\n".format(options_txt_path))
        print("{}\n".format(options_string))

    return options_string

def error(msg):
    sys.stderr.write("ERROR: {}\n".format(msg))
    sys.stderr.flush()
    sys.exit(1)

def execute(command):

    #TODO: Figure out some way to replicate this with plumbum, since plumbum is easier?
    print('Executing {}'.format(command) + os.linesep)
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               shell=True)
    (stdout, stderr) = process.communicate()
    print(stdout.decode() + os.linesep)
    print(stderr.decode() + os.linesep)

#############################
# Script-specific Functions #
#############################

def parse_metadata(metadata_file):

    if os.path.isfile(metadata_file):
        metadata_df = pd.read_table(metadata_file,delimiter='\t',header=0,comment='#')
    else:
        error("Metadata file {} can not be found".format(metadata_file))
    
    if args.debug:
        print("These are the column headings for {}:\n".format(metadata_file))
        print(list(metadata_df))

    return metadata_df

def run_centrifuge(reads, options):
    Status = ''

    #do stuff

    return Status

def run_rna_align(reads, options):
    Status = ''

    # do stuff

    return Status

def run_htseq(alignments, options):
    Status = ''

    # do stuff

    return Status

#main

