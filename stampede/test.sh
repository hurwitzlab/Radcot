#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 68
#SBATCH -t 2:00:00
#SBATCH -p development
#SBATCH -J radcot
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user scottdaniel@email.arizona.edu

export INPUT_DIR="$WORK/radcot-data"
#export OUT_DIR="$WORK/quick-test"
#export GENOME_DIR="$WORK/gtemp"
export OUT_DIR="$WORK/test-of-radcot"
export GENOME_DIR="$WORK/genomes"

bash run.sh -o $OUT_DIR \
    -i $INPUT_DIR \
    -g $GENOME_DIR \
    -m $WORK/subset_crc_mouse_metadata_simreps.txt \
    --debug --threads 12 \
    --procs 4 \
    --skip-centrifuge \
    --skip-rna-align

