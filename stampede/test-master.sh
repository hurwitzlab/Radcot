#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -p development
#SBATCH -J cntrfuge
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user scottdaniel@email.arizona.edu

OUT_DIR="$WORK/centrifuge_out"

export MY_PARAMRUN="$HOME/launcher/paramrun"

[[ -d "$OUT_DIR" ]] && rm -rf $OUT_DIR/*

bash run.sh -q "$WORK/in/dna" -f fastq -o $OUT_DIR -x 9606,32630
#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 12
#SBATCH -t 00:30:00
#SBATCH -p development
#SBATCH -J bowtie2
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user scottdaniel@email.arizona.edu

#for local testing#####
#if the singularity.conf is right, then /vagrant should be auto-shared
export WORK="/vagrant"
#export GUEST="/work"
########################

export OUT_DIR="$WORK/cancer_out"

#export MY_PARAMRUN="$HOME/launcher/paramrun"

[[ -d "$OUT_DIR" ]] && rm -rf $OUT_DIR/*

#    -g "$WORK/genomes" \

bash run.sh \
    -x "$WORK/bt2_index/genome" \
    -1 "$WORK/rna/cancer/RNA_cancer_R1_sample_01.fastq.gz $WORK/rna/cancer/RNA_cancer_R1_sample_02.fastq.gz $WORK/rna/cancer/RNA_cancer_R1_sample_03.fastq.gz" \
    -2 "$WORK/rna/cancer/RNA_cancer_R2_sample_01.fastq.gz $WORK/rna/cancer/RNA_cancer_R2_sample_02.fastq.gz $WORK/rna/cancer/RNA_cancer_R2_sample_03.fastq.gz" \
    -O $OUT_DIR \
    -f fastq -t 4 \
    -a end-to-end \
    -e sensitive
#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -p development
#SBATCH -J cntrfuge
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user scottdaniel@email.arizona.edu

#for local testing#####
#if the singularity.conf is right, then /vagrant should be auto-shared
#export WORK="/vagrant"
export WORK="$HOME/singularity-vm"
export GFF_DIR="$wORK/gffs" #also has genome.fa
########################

export OUT_DIR="$WORK/cuffdiff_test"

#export MY_PARAMRUN="$HOME/launcher/paramrun"

[[ -d "$OUT_DIR" ]] && rm -rf $OUT_DIR/*

#-i "$WORK/genomes"

cuffdiff $GFF_DIR/transcripts.gtf \
    -C sample-sheet.txt \
    -o $OUT_DIR \
    -p 4
