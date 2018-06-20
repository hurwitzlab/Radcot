#!/bin/bash

#SBATCH -J radcot
#SBATCH -A iPlant-Collabs 
#SBATCH -N 1
#SBATCH -n 68
#SBATCH -t 24:00:00
#SBATCH -p normal
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user scottdaniel@email.arizona.edu

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

#check for unset (i.e. [blank]) variables
set -u

export MAINIMG="radcot.img"

if [[ ! -e $MAINIMG ]]; then
    echo "Need the singularity image to work!"
    echo "Go into the /singularity dirs and \"make img\"!"
    exit 1
fi

function HELP() {
    singularity run $MAINIMG -h
    exit 0
}

#
# Show HELP if no arguments
#
[[ $# -eq 0 ]] && echo "Need some arguments" && HELP

#else run the MASTER script
singularity run $MAINIMG $@

echo "Done, look in OUT_DIR \"$OUT_DIR\""
echo "Comments to Scott Daniel <scottdaniel@email.arizona.edu>"

