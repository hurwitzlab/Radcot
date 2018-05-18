#!/usr/bin/env python3

###############################################################################
#                                                                             #
#    cuffKeggR.py                                                             #
#                                                                             #
#    A wrapper script that uses cuffdiff to compare RNA counts between        #
#    samples at a species level, generates graphs with R                      #
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

import glob
import os
import sys
import subprocess
import argparse
from pprint import pprint

#WORK env var will be present on TACC
#But may not be set when testing locally
if os.getenv('WORK') is None:
    os.environ['WORK'] = './'

####################
# ARGUMENTS ########
####################

parser = argparse.ArgumentParser(description=
        "The script essentially wraps cuffdiff for comparing\n"
        "RNA counts between samples.\n",
        formatter_class=argparse.RawTextHelpFormatter)

inputs = parser.add_argument_group('Required Inputs and Parameters',
        "<m1>, <m2>, <r> can be comma-separated lists (no whitespace)\n"
        "and can be specified many times.\n"
        "E.g. -U file1.fq,file2.fq -U file3.fq.")

inputs.add_argument('-i', '--input-dir', 
        dest='input_dir', metavar='DIRECTORY', 
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

gen_opts.add_argument('-k', '--keep-sam', dest='keep_sam', 
        action='store_true',
        help="If enabled, SAM files will be preserved \n"
        "during BAM file generation. Without this option,\n"
        "THERE WILL BE NO SAM FILES.")

gen_opts.add_argument('-n', '--sam-name', 
        dest='sam_name', metavar='FILENAME', 
        default='bowtie2-run.sam',
        help="Filename to use for output sam/bam.")

gen_opts.add_argument('-l', '--log-file', dest='log_fn', 
        metavar='FILENAME', default='bowtie2-read-mapping.log',
        help="Log file name")

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

bowtie2_opts.add_argument('--additional',
        dest='more_args', metavar='STRING',
        default='',
        help="Additional arguments to pass to bowtie2.\n"
        "Put options in single quotes to protect from shell.")

args = parser.parse_args()

####################
# SETS AND CHECKS ##
####################

#check for args that need to be set (the app.json / agave api should do this too)

###############
# FUNCTIONS ###
###############

def error(msg):
    sys.stderr.write("ERROR: {}\n".format(msg))
    sys.stderr.flush()
    sys.exit(1)

def cat_gff(gff_dir,gff):
    file_list = glob.glob(gff_dir + "/*.gff")
    with open(gff+'.gff','w') as w_file:
        for filen in file_list:
            with open(filen, 'rU') as o_file:
                write(o_file, w_file)

        return w_file.name


def execute(command, logfile):

    logfile.write('Executing {}'.format(command) + os.linesep)
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               shell=True)
    (stdout, stderr) = process.communicate()

    logfile.write(stderr + os.linesep)

def make_graphs():

    return None

def cuffdiff(gff):
    
    inFmt = False
    if 'fasta' == args.input_fmt:
        inFmt = '-f'
    if 'fastq' == args.input_fmt:
        inFmt = '-q'
    if 'fq' == args.input_fmt:
        inFmt = '-q'

    if not inFmt:
        error('No format selected during bowtie2 search.')

    #i wish python had case switching! ... this is probably more complicated than it needs to be!
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

    sam_out = os.path.join(args.out_dir, args.sam_name)

    return [(bowtie2_cmd, sam_out)]

def to_sam(cmd2run, keep_sam, logfile):

    processCall = ''

    if len(cmd2run) > 1:  # Alternatively, if sam_output

        for (bowtie2, sam_out) in cmd2run:

            bam_out = sam_out.replace('sam', 'bam')

            # Keep sam file or go directly to bam. Two ways of keeping sam file... one of which sometimes hasn't worked

            if keep_sam:
                processCall = '{1} -S {2} && samtools view --threads {0} -Sb -o {3} {2}'.format(args.threads, bowtie2, sam_out, bam_out)

            else:
                processCall = '{1} | samtools view --threads {0} -Sb - > {2}'.format(args.threads, bowtie2, bam_out)

            execute(processCall, logfile)

    if len(cmd2run) == 1:  # just double-checking

        bowtie2, sam_out = cmd2run[0]
        bam_out = sam_out.replace('sam', 'bam')

        if keep_sam:
            processCall = '{1} -S {2} && samtools view --threads {0} -Sb -o {3} {2}'.format(args.threads, bowtie2, sam_out, bam_out)

        else:
            processCall = '{1} | samtools view --threads {0} -Sb - > {2}'.format(args.threads, bowtie2, bam_out)
            #output = subprocess.check_output(processCall, shell=True)  # Trusted input only!

        execute(processCall, logfile)  # bowtie uses stderr to print out args

    if (len(cmd2run) < 1) and (len(cmd2run) != 1):  # This *should not* ever be triggered
        logfile.write('ERROR: Not >1 or ==1, but sam file specified')

##################
# THE MAIN LOOP ##
##################

if __name__ == '__main__':
   
    #make the out dir if does not exist
    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    #make the log file
    args.log_fn = os.path.join(args.out_dir,args.log_fn)
    log = open(args.log_fn, 'w', 0)  # Set buffer size to 0 to force flushing to disk

    #DEBUG#
    if args.debug:
        log.write('all the arguments:' + os.linesep)
        pprint(args, log)

        log.write('.gff files in {}:'.format(args.gff_dir) + os.linesep)
        pprint(glob.glob(args.gff_dir,"*.gff"),log)
    #END DEBUG#
    
    if os.path.isfile(args.bt2_idx + '.1.bt2') or os.path.isfile(args.bt2_idx + '.1.bt2l'):
        log.write('Bowtie2 index, {}, already exists... assuming its ok'.format(args.bt2_idx) + os.linesep)
        bt2_db_base = args.bt2_idx
    else:
        bt2_db_base = prepare_bowtie_db(args.input_dir, args.bt2_idx, log)

    log.write('Bowtie2 base db: {}'.format(bt2_db_base) + os.linesep)

    cmd_and_sam = bowtie(bt2_db_base)

    to_sam(cmd_and_sam, args.keep_sam, log)

    log.write('Program Complete, Hopefully it Worked!')
    log.close()
