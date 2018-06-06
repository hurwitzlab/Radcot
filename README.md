# Radcot
Determines bacterial genomes from DNA reads then shows gene expression differences from RNA reads

This is a monolithic program that combines the three steps
(using the programs centrifuge, bowtie2, samtools, htseq-count, deseq2, and numerous python and R libraries)

### [centrifuge-patric - Radcot Part One - 1st step](https://github.com/hurwitzlab/centrifuge-patric)
- Identify bacterial species from a metagenomic sample with DNA alignment
- Download genomes and annotations of said species

### [bowtie-samtools - Radcot Part Two - 2nd step](https://github.com/hurwitzlab/bowtie-samtools)
- Align RNA reads to abundanct bacterial genomes

### [count-deseq - Radcot Part Three - 3rd step](https://github.com/hurwitzlab/count-deseq)
- Get counts of RNA alignments to transcript-producing genes (CDS)
- Run Deseq2 guided by the experiment setup file (metadata.txt) that will output graphs / tables of differentially expressed genes

## How to use:
1. Use https://www.imicrobe.us/#/apps to access the app with your [cyverse login](http://www.cyverse.org/create-account)
OR
1. `git clone https://github.com/hurwitzlab/radcot` on a linux system where you have admin rights
2. `git pull && git submodule update --init --recursive` #Gets the individual steps
3. Install [singularity](http://singularity.lbl.gov/all-releases)
4. Build the singularity image in each `/singularity/` by executing `make img`
5. Upload all the files to an HPC with a slurm scheduler<sup>1</sup>
6. Run the program with `sbatch run.sh [arguments]`

---
<sup>1</sup>I assume this can be adapted to run on other 
batch-scheduled high-performance computer systems 
but this has not been tested.
