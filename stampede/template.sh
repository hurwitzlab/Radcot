#!/bin/bash

echo "in_dir \"${in_dir}\""
echo "metadata \"${metadata}\""
echo "cent_opts \"${cent_opts}\""
echo "patric_opts \"${patric_opts}\""
echo "bowtie2_opts \"${bowtie2_opts}\""
echo "htseq_count_opts \"${htseq_count_opts}\""
echo "deseq2_opts \"${deseq2_opts}\""
echo "out_dir \"${out_dir}\""
echo "genome_dir \"${genome_dir}\""
echo "bt2_idx \"${bt2_idx}\""
echo "debug \"${debug}\""
echo "threads\"${threads}\""
echo "procs\"${procs}\""
echo "skip_cent \"${skip_cent}\""
echo "skip_rna \"${skip_rna}\""

bash run.sh ${in_dir} ${metadata} ${cent_opts} ${patric_opts} ${bowtie2_opts} \
    ${htseq_count_opts} ${deseq2_opts} ${out_dir} ${genome_dir} ${bt2_idx} \
    ${debug} ${threads} ${procs} ${skip_cent} ${skip_rna} 
