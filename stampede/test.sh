#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 68
#SBATCH -t 06:00:00
#SBATCH -p normal
#SBATCH -J radcot
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user scottdaniel@email.arizona.edu

#for local testing#####
#if the singularity.conf is right, then /vagrant should be auto-shared
if [[ -z $MODULEPATH ]]; then
    export WORK="/vagrant"
fi

export INPUT_DIR="$WORK/in-smaller"
########################

export OUT_DIR="$WORK/radcot_test"
export GENOME_DIR="$WORK/genomes"

#export MY_PARAMRUN="$HOME/launcher/paramrun"

#[[ -d "$OUT_DIR" ]] && rm -rf $OUT_DIR/*

#-i "$WORK/genomes"

bash run.sh -o $OUT_DIR \
    -i $INPUT_DIR \
    -g $GENOME_DIR \
    -m $WORK/subset_crc_mouse_metadata_simreps.txt \
    --debug --threads 12 \
    --procs 4 \
    --skip-centrifuge \
    --skip-rna-align

#from makefile in ../scripts:
#test_runall:
#	singularity exec $(IMG) runAll.py \
#		-i $(WORK)/in-smaller \
#		-o $(WORK)/radcot-test \
#		-m $(WORK)/subset_crc_mouse_metadata.txt \
#		-C $(WORK)/centrifuge-opts.txt \
#		-p $(WORK)/patric-opts.txt \
#		-B $(WORK)/bowtie2-opts.txt \
#		-H $(WORK)/htseq-opts.txt \
#		-D $(WORK)/deseq2-opts.txt \
#		-d -t 4 -P 4
