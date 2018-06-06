#!/bin/bash

#SBATCH -J radcot
#SBATCH -A iPlant-Collabs 
#SBATCH -N 12
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -p normal

# Author: Scott G. Daniel <scottdaniel@email.arizona.edu>

###Uncomment when back on tacc#
#echo "#### Current modules after app.json processing:"
#module list 2>&1
#echo "#### LOADING TACC-SINGULARITY ####"
#module load tacc-singularity 2>&1
#echo "#### LOADING LAUNCHER ####"
#module load launcher 2>&1
#echo "#### Current modules after run.sh processing:"
#module list 2>&1
#
# Set up defaults for inputs, constants
#

#check for unset (i.e. [blank]) variables
set -u

#If this repo is properly checked out, this will work#
#Otherwise, you probably forgot to: #
#git pull && git submodule update --init --recursive#

STEPONE="../../01-centrifuge-patric"
STEPTWO="../../02-bowtie-samtools"
STEPTHREE="../../03-count-deseq"

if [[ ! -d $STEPONE || ! -d $STEPTWO || ! -d $STEPTHREE ]]; then
    echo "Can not find the required submodules"
    echo "You need to \"git pull && git submodule update --init --recursive\""
    exit 1
fi

CENTIMG="$STEPONE/stampede/centrifuge-patric.img"
BOWTIMG="$STEPTWO/stampede/bowtie-samtools.img"
HTSQIMG="$STEPTHREE/stampede/count-deseq.img"

OUT_DIR="$PWD/radcot-out"
#
# Some needed functions
#
function lc() { 
    wc -l "$1" | cut -d ' ' -f 1 
}

function HELP() {
    echo "https://github.com/hurwitzlab/Radcot/tree/master"
    exit 0
}

#
# Show HELP if no arguments
#
[[ $# -eq 0 ]] && echo "Need some arguments" && HELP

echo "Done, look in OUT_DIR \"$OUT_DIR\""
echo "Comments to Scott Daniel <scottdaniel@email.arizona.edu>"

