#!/usr/bin/env python

###############################################################################
#                                                                             #
#    patric_bowtie2.py                                                        #
#                                                                             #
#    A wrapper script for Bowtie2, that runs Bowtie and Samtools on           #
#    read files using a directory of PATRIC (patricbrc.org) genomes           #
#                                                                             #
#    Copyright (C) Benjamin Bolduc, Scott G Daniel                            #
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

__author__ = ["Ben Bolduc", "Scott G Daniel"]
__copyright__ = "Copyright 2018"
__credits__ = ["Ben Bolduc","Scott G Daniel"]
__license__ = "LGPLv3"
__maintainer__ = ["Ben Bolduc", "Scott G Daniel"]
__email__ = ["bbolduc.chem@gmail.com","scottdaniel@email.arizona.edu"]
__status__ = "Development"

import glob
import os
import sys
import subprocess
import argparse
import itertools
from pprint import pprint
from Bio import SeqIO

#WORK env var will be present on TACC
#But may not be set when testing locally
if os.getenv('WORK') is None:
    os.environ['WORK'] = './'

####################
# ARGUMENTS ########
####################

parser = argparse.ArgumentParser(description=
        "The script essentially wraps bowtie2 for aligning\n"
        "reads against the same reference sequence collection.\n",
        formatter_class=argparse.RawTextHelpFormatter)

inputs = parser.add_argument_group('Required Inputs and Parameters',
        "<m1>, <m2>, <r> can be comma-separated lists (no whitespace)\n"
        "and can be specified many times.\n"
        "E.g. -U file1.fq,file2.fq -U file3.fq.")

inputs.add_argument('-g', '--genome-dir', 
        dest='genome_dir', metavar='DIRECTORY', 
        default=os.path.join(os.getenv('WORK'),'genomes'),
        help="The Directory containing individual genomes\n"
        "that will be pasted together. The created genome.fna\n"
        "will be indexed for bowtie2")

inputs.add_argument('-x', '--bt2-idx', 
        dest='bt2_idx', metavar='FILENAME', 
        default=os.path.join(os.getenv('WORK'),'bt2_index/','genome'),
        help="Index filename prefix (minus trailing .X.bt2).\n"
        "This will also be the name of the fasta file, E.g. [bt2-idx].fna\n"
        "NOTE: Bowtie 1 and Bowtie 2 indexes are not compatible.")

inputs.add_argument('-1', '--m1', 
        dest='reads_forward', metavar='STRING',
        default='',
        help="Files with #1 mates, paired with files in <m2>.\n"
        "Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).\n")

inputs.add_argument('-2', '--m2', 
        dest='reads_reverse', metavar='STRING',
        default='',
        help="Files with #2 mates, paired with files in <m1>.\n"
        "Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).\n")

inputs.add_argument('-U', '--unpaired', 
        dest='reads_unpaired', metavar='STRING',
        default='',
        help="Files with unpaired reads.\n"
        "Could be gzip'ed (extension: .gz) or bzip2'ed (extension: .bz2).")

inputs.add_argument('-f', '--fmt', '--input-format', 
        dest='input_fmt', 
        choices=['fastq', 'fasta', 'fq'],
        default='fastq',
        help="File format of reads to be aligned. \n"
        "Compressed files (*.gz, *.tar.gz) will be automatically\n"
        "recognized.")

gen_opts = parser.add_argument_group('General Options')  

gen_opts.add_argument('-O', '--out-dir', dest='out_dir', type=str,
        default=os.getcwd, help="Output directory to put all the\n"
        "results in.")

gen_opts.add_argument('-n', '--bam-name', 
        dest='bam_name', metavar='FILENAME', 
        default='bowtie2-run.bam',
        help="Filename to use for output bam. \n"
        "This is usually defined by the calling script \n"
        "based on the input read names")

gen_opts.add_argument('-l', '--log-file', dest='log_fn', 
        metavar='FILENAME', default='bowtie2-read-mapping.log',
        help="Log file name")
###
bowtie2_opts = parser.add_argument_group('Bowtie2 Alignment Options')

bowtie2_opts.add_argument('-a', '--alignment-type', 
        dest='align_type', choices=['end-to-end', 'local'], 
        default='end-to-end',
        help="Whether the entire read must align \n"
        "(end-to-end) or only a local region (local).")

bowtie2_opts.add_argument('-e', '--end-to-end-presets', 
        dest='global_presets', metavar='STRING',
        choices=['very-fast', 'fast', 'sensitive', 'very-sensitive'],
        default="sensitive", 
        help="Presets for end-to-end alignments:\n"
        "very-fast: -D 5 -R 1 -N 0 -L 22 -i S,0,2.50\n"
        "fast: -D 10 -R 2 -N 0 -L 22 -i S,0,2.50\n"
        "sensitive: -D 15 -R 2 -N 0 -L 22 -i S,1,1.15 (default)\n"
        "very-sensitive: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50")

bowtie2_opts.add_argument('-L', '--local-presets', 
        dest='local_presets', metavar='STRING',
        choices=['very-fast-local', 'fast-local', 
            'sensitive-local', 'very-sensitive-local'],
        default='sensitive-local',
        help="Presets for local alignments:\n"
        "very-fast-local: -D 5 -R 1 -N 0 -L 25 -i S,1,2.00\n"
        "fast-local: -D 10 -R 2 -N 0 -L 22 -i S,1,1.75\n"
        "sensitive-local: -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)\n"
        "very-sensitive-local: -D 20 -R 3 -N 0 -L 20 -i S,1,0.50")

bowtie2_opts.add_argument('-N', '--non-deterministic', 
        dest='non_deterministic', action='store_true',
        help="Bowtie 2 will use the current time to \n"
        "re-initialize the pseudo-random number generator.\n"
        "Useful when the input consists of many identical reads.")

bowtie2_opts.add_argument('-5', '--trim5', 
        dest='trim5', metavar='INT', 
        type=int, default=0,
        help="Trim X bases from 5'/left end of reads.")

bowtie2_opts.add_argument('-3', '--trim3', 
        dest='trim3', metavar='INT', 
        type=int, default=0,
        help="Trim X bases from 3'/right end of reads.")

bowtie2_opts.add_argument('-I', '--minins', 
        dest='minins', metavar='INT', 
        type=int, default=0,
        help="minimum fragment length (0)")

bowtie2_opts.add_argument('-X', '--maxins', 
        dest='maxins', metavar='INT', 
        type=int, default=2000,
        help="maximum fragment length (2000)")

bowtie2_opts.add_argument('-t', '--threads', 
        dest='threads', metavar='INT', 
        type=int, default=1,
        help="number of alignment threads to launch (1)")

bowtie2_opts.add_argument('-A', '--additional',
        dest='more_args', metavar='STRING',
        default='',
        help="Additional arguments to pass to bowtie2.\n"
        "Put options in single quotes to protect from shell.")

args = parser.parse_args()

####################
# SETS AND CHECKS ##
####################


#check for args that need to be set (the app.json / agave api should do this too)
preset = False
if args.align_type == "end-to-end":
    preset = args.global_presets
if args.align_type == "local":
    preset = args.local_presets

if args.input_fmt not in ['fasta', 'fastq', 'fq']:
    error('ERROR: Input format must be either fasta or fastq formatted')

###############
# FUNCTIONS ###
###############

def error(msg):
    sys.stderr.write("ERROR: {}\n".format(msg))
    sys.stderr.flush()
    sys.exit(1)

def cat_fasta(genome_dir,bt2_idx):
    file_list = glob.glob(genome_dir + "/*.fna")
    with open(bt2_idx+'.fna','w') as w_file:
        for filen in file_list:
            with open(filen, 'rU') as o_file:
                seq_records = SeqIO.parse(o_file, 'fasta')
                SeqIO.write(seq_records, w_file, 'fasta')

        return w_file.name


def execute(command, logfile):

    logfile.write('Executing {}'.format(command) + os.linesep)
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               shell=True)
    (stdout, stderr) = process.communicate()

    logfile.write(stderr + os.linesep)

def prepare_bowtie_db(genome_dir, bt2_idx, logfile):

    bt2_db_fasta = args.bt2_idx + '.fna'
    db_dir = os.path.dirname(args.bt2_idx) #E.g. /vagrant/bt2_idx
    # Ensure there's a genome.fna file or make one if not
    if not os.path.isfile(bt2_db_fasta): #E.g. /vagrant/bt2_idx/genome.fna
        pprint("No input bowtie2 db specified," \
               + " concatenating fastas in {}".format(args.genome_dir),logfile)
        #make the directory for the bt2 index if not there
        if not os.path.isdir(db_dir): 
            os.makedirs(db_dir) #E.g. mkdir /vagrant/bt2_idx
        bt2_db_fasta = cat_fasta(args.genome_dir,args.bt2_idx)
        pprint("Created a combined genome for you: {}".format(bt2_db_fasta),logfile)

    bowtie_db_cmd = 'bowtie2-build --threads {} -f {} {}'.format(args.threads, bt2_db_fasta, args.bt2_idx)
    execute(bowtie_db_cmd, logfile)

    return args.bt2_idx


def bowtie(bowtie2_db):
    
    inFmt = False
    if 'fasta' == args.input_fmt:
        inFmt = '-f'
    if 'fastq' == args.input_fmt:
        inFmt = '-q'
    if 'fq' == args.input_fmt:
        inFmt = '-q'

    if not inFmt:
        error('No format selected during bowtie2 search.')

    #this logic should already have been done in the run.sh
    #but i guess we are double-checking?
    if args.reads_unpaired and args.reads_forward and args.reads_reverse:
        input_cmd = '-1 {} -2 {} -U {}'.format(args.reads_forward,args.reads_reverse,args.reads_unpaired)
    elif args.reads_unpaired and not (args.reads_forward and args.reads_reverse):
        input_cmd = '-U {}'.format(args.reads_unpaired)
    elif (args.reads_forward and args.reads_reverse) and not args.reads_unpaired:
        input_cmd = '-1 {} -2 {}'.format(args.reads_forward,args.reads_reverse)
    else:
        error("Something is wrong in how the input string is formatted\n"
        "Did you only enter forward reads? Did you only enter reverse reads?\n")

    bowtie2_cmd = 'bowtie2 {} --phred33 --{} --{} -p {} -I {} -X {} --no-unal {} -x {} {}'.format(
        inFmt, args.align_type, preset, args.threads, args.minins, 
        args.maxins, args.more_args, bowtie2_db, input_cmd)

    bam_out = os.path.join(args.out_dir, args.bam_name)

    return [(bowtie2_cmd, bam_out)]

def to_bam(cmd2run, logfile):

    processCall = ''

    for (bowtie2, bam_out) in cmd2run:
        
        #DEBUG#
        log.write('Running bowtie2 and converting to bam' + os.linesep)

        convert_to_bam = '{} | samtools view --threads {} -bT {} - > {}'.format( bowtie2, args.threads, args.bt2_idx + '.fna', bam_out + '.tmp')

        execute(convert_to_bam, logfile)

        #Debug#
        log.write('Sorting bam by position' + os.linesep)
        
        sort_bam = 'samtools sort --threads {} -n {} > {}'.format(args.threads, bam_out + '.tmp', bam_out)

        execute(sort_bam, logfile)
        
        if os.path.isfile(bam_out):
            os.remove(bam_out + '.tmp')

##################
# THE MAIN LOOP ##
##################

if __name__ == '__main__':
   
    #make the out dir if does not exist
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    #make the log file
    args.log_fn = os.path.join(args.out_dir,args.log_fn)
    log = open(args.log_fn, 'w', 0)  # Set buffer size to 0 to force flushing to disk

    #DEBUG#
#    log.write('ALL THE ARGUMENTS:' + os.linesep)
#    pprint(args, log)
#
#    log.write('Directory contents for genomes:' + os.linesep)
#    pprint(os.listdir(args.genome_dir),log)
    #END DEBUG#
    
    if os.path.isfile(args.bt2_idx + '.1.bt2') or os.path.isfile(args.bt2_idx + '.1.bt2l'):
        log.write('Bowtie2 index, {}, already exists... assuming its ok'.format(args.bt2_idx) + os.linesep)
        bt2_db_base = args.bt2_idx
    else:
        bt2_db_base = prepare_bowtie_db(args.genome_dir, args.bt2_idx, log)

    log.write('Bowtie2 base db: {}'.format(bt2_db_base) + os.linesep)

    cmd_and_bam = bowtie(bt2_db_base)

    to_bam(cmd_and_bam, log)

    log.write('Program Complete, Hopefully it Worked!')
    log.close()
