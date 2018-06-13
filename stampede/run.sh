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
if [[ $MODULEPATH ]]; then
    echo "#### LOADING TACC-SINGULARITY ####"
    module load tacc-singularity 2>&1
    echo "#### LOADING LAUNCHER ####"
    module load launcher 2>&1
else
    echo "No modulepath"
fi
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

export PRJROOT=".."
STEPONE="$PRJROOT/01-centrifuge-patric"
STEPTWO="$PRJROOT/02-bowtie-samtools"
STEPTHREE="$PRJROOT/03-count-deseq"

if [[ ! -d $STEPONE || ! -d $STEPTWO || ! -d $STEPTHREE ]]; then
    echo "Can not find the required submodules"
    echo "You need to \"git pull && git submodule update --init --recursive\""
    exit 1
fi

export MAINIMG="radcot.img"
export CENTIMG="$STEPONE/stampede/centrifuge-patric.img"
export BOWTIMG="$STEPTWO/stampede/bowtie-sam.img"
export HTSQIMG="$STEPTHREE/stampede/count-deseq.img"

if [[ ! -e $MAINIMG || ! -e $CENTIMG || ! -e $BOWTIMG || ! -e $HTSQIMG ]]; then
    echo "Need the singularity images to work!"
    echo "Go into the /singularity dirs and \"make img\"!"
    exit 1
fi

OUT_DIR="$PWD/radcot-out"
#
# Some needed functions
#
function lc() { 
    wc -l "$1" | cut -d ' ' -f 1 
}

function HELP() {
    singularity run radcot.img -h
    exit 0
}

#
# Show HELP if no arguments
#
[[ $# -eq 0 ]] && echo "Need some arguments" && HELP

#else run the MASTER script
singularity run radcot.img $@

echo "Done, look in OUT_DIR \"$OUT_DIR\""
echo "Comments to Scott Daniel <scottdaniel@email.arizona.edu>"

