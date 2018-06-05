#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 68
#SBATCH -t 02:00:00
#SBATCH -p development
#SBATCH -J deseq2
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user scottdaniel@email.arizona.edu

#echo "Current working directory is $(pwd)"
#export CWD=$(pwd)
#for local testing#####
#if the singularity.conf is right, then /vagrant should be auto-shared
#export WORK="/vagrant"
export BAMS_DIR="$WORK/bams"
export GFF_DIR="$WORK/genomes" #also has genome.fa
export METADATA="metadata.txt"
########################


export OUT_DIR="$WORK/deseq_test"
[[ ! -d "$OUT_DIR" ]] && mkdir -p "$OUT_DIR"
#export MY_PARAMRUN="$HOME/launcher/paramrun"

#[[ -d "$OUT_DIR" ]] && rm -rf $OUT_DIR/*

#-i "$WORK/genomes"
#
./count-deseq.py -G $GFF_DIR \
    -g $BAMS_DIR/all.RefSeq.gff \
    -d -t 68 \
    --skip-counting \
    -b $BAMS_DIR \
    -m $METADATA \
    -o $OUT_DIR \
    -C htseq_count_options.txt \
    -D deseq2_options.txt
