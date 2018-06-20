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

#important env vars
#CENTIMG="centrifuge-patric.img"
#BOWTIMG="bowtie-sam.img"
#HTSQIMG="count-deseq.img"

import os, sys, argparse, glob, subprocess
import pandas as pd
from pprint import pprint

#WORK env var will be present on TACC
#But may not be set when testing locally
if os.getenv('WORK') is None:
    os.environ['WORK'] = './' #this is how we could set launcher variables, etc.

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

inputs.add_argument('-i', '--in-dir',
        dest='in_dir', metavar='DIRECTORY',
        default=os.path.join(os.getenv('WORK'),'in'),
        help="Input directory with all the dna / rna reads.\n"
        "This will be prepended to each file specified in\n"
        "the metadata file, so you don't need to specifiy it\n"
        "there too. [ Default = $WORK/in ]")

inputs.add_argument('-o', '--out-dir', 
        dest='out_dir', metavar='DIRECTORY',
        default=os.getcwd(),
        help="Output directory to put all the\n"
        "results in. [ Default = Current working dir ]")

inputs.add_argument('-m', '--metadata', 
        dest='metadata', metavar='FILENAME',
        default=os.path.join(os.getenv('WORK'),'metadata.txt'),
        help="File containing file / sample information.\n"
        "Use the metadata_template to start and DO NOT change\n"
        "headings. [ Default = $WORK/metadata.txt ]")

step_one = parser.add_argument_group('Step one option\'s files (01-centrifuge-patric)')

step_one.add_argument('-C', '--centrifuge-opts', 
        dest='cent_opts', metavar='FILENAME', 
        default=os.path.join(os.getenv('WORK'),'centrifuge-opts.txt'),
        help="File with additional options for centrifuge\n"
        "Please format your options like so:\n"
        "-o1 option1\n"
        "-o2 option2\n"
        "[ Default = $WORK/centrifuge-opts.txt ]")

step_one.add_argument('-p', '--patric-opts', 
        dest='patric_opts', metavar='FILENAME', 
        default=os.path.join(os.getenv('WORK'),'patric-opts.txt'),
        help="File with additional options for getting genomes from patric\n"
        "Please format your options like so:\n"
        "-o1 option1\n"
        "-o2 option2\n"
        "[ Default = $WORK/patric-opts.txt ]")

step_two = parser.add_argument_group(
        'Step two option\'s files (02-bowtie-samtools)')

step_two.add_argument('-B', '--bowtie2-opts', 
        dest='bowtie2_opts', metavar='FILENAME',
        default=os.path.join(os.getenv('WORK'),'bowtie2-opts.txt'),
        help="File with additional options for bowtie2 alignment for RNA\n"
        "Please format your options like so:\n"
        "-o1 option1\n"
        "-o2 option2\n"
        "[ Default = $WORK/bowtie2-opts.txt ]")

step_three = parser.add_argument_group(
        'Step three option\'s files (03-count-deseq)')

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
        help="Number of alignment threads to launch.\n"
        " [ Default = 1 ] ")

gen_opts.add_argument('-P', '--procs',
        dest='procs', metavar='INT',
        type=int, default=1,
        help="Number of parallel processes to launch\n"
        "with gnu-parallel (this multiplies with threaded programs).\n"
        " [ Default = 1 ] ")

gen_opts.add_argument('--skip-centrifuge',
        dest='skip_cent', 
        action='store_true',
        help="Skip step one, centrifuge.\n"
        "NOTE: this means that you already have genomes\n"
        "ready to align your RNA reads to.\n")

gen_opts.add_argument('--skip-rna-align',
        dest='skip_rna',
        action='store_true',
        help="Skip steps one and two.\n"
        "This means you already have genomes / gffs AND\n"
        "sam files that are sorted by name for htseq_count.\n")

args = parser.parse_args()

#######################
# GENERAL FUNCTIONS ###
#######################

#this loads modules on TACC or UA HPC
#also has the side-effect of reloading the entire environment so CAUTION
#currently, this isn't being used because we have the ../stampede/run.sh
#that does the module loading for us
def module_load(module_name):
     #load the module and print out all environmental variables
    command = 'module load ' + module_name + ' && env'
    proc = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)

    #parse through each variable and set it in the working environment of this
    #script, CAUTION: this will overrite anything there
    for aline in proc.stdout:
        (key, _, value) = aline.decode().partition("=")
        os.environ[key] = value.rstrip('\n')
    proc.communicate()

#basic checker, check that options file exists and then check that each line begins with a '-', then parse
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

    print('Executing {}'.format(command) + os.linesep)
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               shell=True)
    (stdout, stderr) = process.communicate()
    print(stdout.decode() + os.linesep)
    print(stderr.decode() + os.linesep)

    return process.returncode

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

def parse_reads(panda_df, panda_column):
    """Simple function that takes a panda column
    and spits out a string that centrifuge and bowtie2 like"""
    read_full_paths = []
    read_string = ''

    for read in list(panda_df[panda_column]):
        if type(read) == float: #if the field is empty, pandas returns a float
            continue
        else:
            read_full_paths.append(os.path.join(args.in_dir,read))
    
    if len(read_full_paths) > 0: #return an empty string if nothing in the column
        read_string = ','.join(read_full_paths)

    return read_string

def run_centrifuge(reads, cent_opts, patric_opts):

    options_string = parse_options_text(cent_opts) + parse_options_text(patric_opts)

    #If we are using this from the metadata txt
    #we will have either paired + unpaired
    #or just paired
    #or just unpaired
    #the -q won't be used
    f_reads = ''
    r_reads = ''
    u_reads = ''

    f_reads = parse_reads(metadata, 'dna_forward')
    r_reads = parse_reads(metadata, 'dna_reverse')
    u_reads = parse_reads(metadata, 'dna_unpaired')

    if args.debug:
        print("These are the reads: forward {}\n".format(f_reads))
        print("reverse {}\n".format(r_reads))
        print("and unpaired {}\n".format(u_reads))

    bin_dir = os.path.dirname(os.path.realpath(__file__))
    cent_script = os.path.join(bin_dir, 'run_centrifuge.py')

    if f_reads and r_reads and not u_reads:

        command = '{} -1 {} -2 {} -o {} {}'.format(cent_script, f_reads,
                    r_reads, args.out_dir, options_string)

        returncode = execute(command)

        if returncode == 0:
            print('{} ran sucessfully, continuing...'.format(cent_script))
        else:
            error('{} failed, exiting'.format(cent_script))

    elif f_reads and r_reads and u_reads:
        
        command = '{} -1 {} -2 {} -U -o {} {}'.format(cent_script, f_reads,
                    r_reads, u_reads, args.out_dir, options_string)

        returncode = execute(command)

        if returncode == 0:
            print('{} ran sucessfully, continuing...'.format(cent_script))
        else:
            error('{} failed, exiting'.format(cent_script))

    elif u_reads and not f_reads or r_reads:

        command = '{} -U {} -o {} {}'.format(cent_script, u_reads,
                    args.out_dir, options_string)

        returncode = execute(command)

        if returncode == 0:
            print('{} ran sucessfully, continuing...'.format(cent_script))
        else:
            error('{} failed, exiting'.format(cent_script))

    else:
        error('No Reads!')
    
def run_rna_align(reads, options):
    Status = ''

    # do stuff

    return Status

def run_htseq(alignments, options):
    Status = ''

    # do stuff

    return Status

##################
# THE MAIN LOOP ##
##################

if __name__ == '__main__':
   
    #make the out dir if does not exist
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    #DEBUG#
    if args.debug:
        print('all the arguments:' + os.linesep)
        pprint(args); print()
    
    #parse metadata
    metadata = parse_metadata(args.metadata)

    #DEBUG: parse all options
#    if args.debug:
#        for key in filter(lambda key: key.endswith('opts'), vars(args)):
#            parse_options_text(vars(args)[key])
    
    #Run centrifuge
    if not args.skip_cent:
        print("Running centrifuge")
        run_centrifuge(metadata, args.cent_opts, args.patric_opts)

    #Run bowtie2
    if not args.skip_rna:
        print("Running bowtie2 alignment for RNA reads")
        run_rna_align(metadata, args.bowtie2_opts)
    
    #Run htseq-count and deseq2
    run_htseq(metadata, args.htseq_count_opts, args.deseq2_opts)

