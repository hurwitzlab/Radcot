# Radcot
Determines bacterial genomes from DNA reads then shows gene expression differences from RNA reads

This is a monolithic program that combines the three steps
(using the programs centrifuge, bowtie2, samtools, htseq-count, deseq2, and numerous python and R libraries)

The idea is that this one program combines all three steps so you do not have to run them separately.

### [centrifuge-patric - Radcot Part One - 1st step](https://github.com/hurwitzlab/centrifuge-patric)
- Identify bacterial species from a metagenomic sample with DNA alignment
- Download genomes and annotations of said species

### [bowtie-samtools - Radcot Part Two - 2nd step](https://github.com/hurwitzlab/bowtie-samtools)
- Align RNA reads to abundanct bacterial genomes

### [count-deseq - Radcot Part Three - 3rd step](https://github.com/hurwitzlab/count-deseq)
- Get counts of RNA alignments to transcript-producing genes (CDS)
- Run Deseq2 guided by the experiment setup file (metadata.txt) that will output graphs / tables of differentially expressed genes

## Installation / Quick start:
1. Use https://www.imicrobe.us/#/apps to access the app with your [cyverse login](http://www.cyverse.org/create-account)
OR
2. `git clone https://github.com/hurwitzlab/radcot` on a linux system where you have admin rights
3. Install [singularity](http://singularity.lbl.gov/all-releases)
4. Build the singularity image in `singularity/` by executing `make img`
5. Upload all the files to an HPC with a slurm scheduler<sup>1</sup>
6. Run the program with `sbatch run.sh [arguments]`

---
<sup>1</sup>This can be adapted to run on other 
batch-scheduled high-performance computer systems 
but this has not been tested.

---

## Example run:
#### Getting the example data
  1. Have [iRODS tools](https://docs.irods.org/master/getting_started/download/) installed on your HPC
  2. Get the example data with `iget -r /iplant/home/shared/imicrobe/example_data/radcot-data`
    * If using iMicrobe, this will be your input directory and you shouldn't need to download
  3. To check you have the data, make sure each DNA and RNA read file matches those found in metadata_template.txt
