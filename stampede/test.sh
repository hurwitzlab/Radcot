#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 02:00:00
#SBATCH -p development
#SBATCH -J deseq2
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user scottdaniel@email.arizona.edu

#for local testing#####
#if the singularity.conf is right, then /vagrant should be auto-shared
export WORK="/vagrant"
export INPUT_DIR="$WORK/in"
########################

export OUT_DIR="$WORK/radcot_test"

#export MY_PARAMRUN="$HOME/launcher/paramrun"

[[ -d "$OUT_DIR" ]] && rm -rf $OUT_DIR/*

#-i "$WORK/genomes"

bash run.sh -o $OUT_DIR \
    -i $INPUT_DIR \
    -m $WORK/subset_crc_mouse_metadata.txt \
    --debug --threads 4
