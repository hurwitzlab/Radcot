#!/bin/bash

#SBATCH -J htseq
#SBATCH -A iPlant-Collabs 
#SBATCH -N 12
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -p normal

# Author: Scott G. Daniel <scottdaniel@email.arizona.edu>

###Uncomment when back on tacc#
#echo "#### Current modules after app.json processing:"
#module list 2>&1
echo "#### LOADING TACC-SINGULARITY ####"
module load tacc-singularity 2>&1
echo "#### LOADING LAUNCHER ####"
module load launcher 2>&1
#echo "#### Current modules after run.sh processing:"
#module list 2>&1
#
# Set up defaults for inputs, constants
#
SING_IMG="count-deseq.img"
OUT_DIR="$PWD/count-deseq-out"
#
# Some needed functions
#
function lc() { 
    wc -l "$1" | cut -d ' ' -f 1 
}

function HELP() {
    singularity exec $SING_IMG count-deseq.py -h
    exit 0
}

#
# Show HELP if no arguments
#
[[ $# -eq 0 ]] && echo "Need some arguments" && HELP

#set -u

#check for centrifuge image
if [[ ! -e "$SING_IMG" ]]; then
    echo "Missing SING_IMG \"$SING_IMG\""
    exit 1
fi

#Run htseq-count and SARTools (which runs deseq2)
singularity exec $SING_IMG count-deseq.py $@

echo "Done, look in OUT_DIR \"$OUT_DIR\""
echo "Comments to Scott Daniel <scottdaniel@email.arizona.edu>"

