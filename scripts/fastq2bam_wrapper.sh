#!/usr/bin/bash

cd /mnt/scratch/steepale/birdman/bermuda/fastq2bam

sample=$1

# Load appropriate modules
module load Python/3.3.2
module load FastQC/0.11.3
module load Trimmomatic/0.33

# Submit job to cluster
qsub -N "fastq2bam_"$sample "python ./scripts/fastq2bam.py $sample"
