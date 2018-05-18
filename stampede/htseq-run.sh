#!/bin/bash

#SBATCH -J bowtie2
#SBATCH -A iPlant-Collabs 
#SBATCH -N 4
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -p normal

# Author: Scott G. Daniel <scottdaniel@email.arizona.edu>

###Uncomment when back on tacc#
#module load tacc-singularity 
#module load launcher
#
# Some needed functions
#
function lc() { 
    wc -l "$1" | cut -d ' ' -f 1 
}

function HELP() {
    
    echo "BAMS="./bowtie2-run.bam" #-b | --bams"
    echo "GFF_DIR="./genomes" #-g | --gff-dir"
    echo "GFF="./genomes/genome.gff" #-G | --gff"
    echo "LOG_FN="cuffdiff.log" #-l | --log-file"
    echo "THREADS="1" #-t | --threads"
    echo "OUT_DIR="./out_dir" #-O | --out-dir"
    echo "RRNA_GFF="./rRNA.gff" #-M | --mask-gff"
    echo "DEBUG="FALSE" #-d | --debug"
    echo "SILENT="FALSE" #-s | --silent"
    echo "MORE_ARGS="" #-A | --additional"

    echo "See cufflinks.py or make_graphs.R for additional help"
    exit 0
}

#
# Show HELP if no arguments
#
[[ $# -eq 0 ]] && echo "Need some arguments" && HELP

set -u

#
# Set up defaults for inputs, constants
#
BAMS="./bowtie2-run.bam" #-b | --bams"
GFF_DIR="./genomes" #-g | --gff-dir"
GFF="./genomes/genome.gff" #-G | --gff"
LOG_FN="cuffdiff.log" #-l | --log-file"
THREADS="1" #-t | --threads"
SING_IMG="cuffKeggR.img" 
OUT_DIR="./out_dir" #-O | --out-dir"
RRNA_GFF="./rRNA.gff" #-M | --mask-gff"
DEBUG="FALSE" #-d | --debug"
SILENT="FALSE" #-s | --silent"
MORE_ARGS="" #-A | --additional"

#Read the arguments
# In case you wanted to check what variables were passed
echo "ARG = $*"

while getopts :b:g:G:l:t:S:O:M:d:s:A:h ARG; do
    case $ARG in
        b)
            BAMS="$OPTARG"
            ;;
        g)
            GFF_DIR="$OPTARG"
            ;;
        G)
            GFF="$OPTARG"
            ;;
        l)
            LOG_FN="$OPTARG"
            ;;
        t)
            THREADS="$OPTARG"
            ;;
        o)
            OUT_DIR="$OPTARG"
            ;;
        M)
            RRNA_GFF="$OPTARG"
            ;;
        d)
            DEBUG=1
            ;;
        s)
            SILENT=1
            ;;
        A)
            MORE_ARGS="$OPTARG"
            ;;
        h)
            HELP
            ;;
        :)
            echo ""$OPTARG" requires an argument"
            ;;
        \?) #unrecognized option - show help
            echo "Invalid option "$OPTARG""
            HELP
            ;;
    esac
done

echo "What it looks like after parsing"
echo "ARG = $*"
#It is common practice to call the shift command 
#at the end of your processing loop to remove 
#options that have already been handled from $@.
shift $((OPTIND -1))


#If you have your own launcher setup on stampede2 just point MY_PARAMRUN at it
#this will override the TACC_LAUNCHER...
#PARAMRUN="${MY_PARAMRUN:-$TACC_LAUNCHER_DIR/paramrun}"

#check for centrifuge image
if [[ ! -e "$SING_IMG" ]]; then
    echo "Missing SING_IMG \"$SING_IMG\""
    exit 1
fi

#
# Verify existence of various directories, files
#
#[[ ! -d "$OUT_DIR" ]] && mkdir -p "$OUT_DIR"
# python script will do these checks

###Uncomment when back on TACC
#if [[ ! -d "$TACC_LAUNCHER_DIR" ]]; then
#    echo "Cannot find TACC_LAUNCHER_DIR \"$TACC_LAUNCHER_DIR\""
#    exit 1
#fi
#
#if [[ ! -f "$PARAMRUN" ]]; then
#    echo "Cannot find PARAMRUN \"$PARAM_RUN\""
#    exit 1
#fi


#Run cuffquant, cuffnorm, and cuffdiff
singularity run $SING_IMG \
-g $GFF_DIR \
-G $GFF \
-l $LOG_FN \
-t $THREADS \
-O $OUT_DIR \
-M $RRNA_GFF \
-d $DEBUG \
-s $SILENT \
-A $MORE_ARGS \
$BAMS

echo "Done, look in OUT_DIR \"$OUT_DIR\""
echo "Comments to Scott Daniel <scottdaniel@email.arizona.edu>"

