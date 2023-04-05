#!/bin/bash

### Trim whole genome sequence reads from 137 Littorin saxatilis samples.
### Reads were squenced by GeneWiz using an Illumina library
### preperation. Adaptor sequences were not clipped.
### This script was designed to run on the Rackham server of UppMax.
### James Reeve - Göteborgs Universitet
### 21/03/2023

### Job parameters
#SBATCH -A snic2022-5-266 
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 4:00:00
#SBATCH --array=1-137%20
#SBATCH -J fastp_trim
#SBATCH --output logfiles/fastp_trim_%a.log 

## Script
module load bioinfo-tools
module load fastp/0.23.2

# File paths
DIR=/proj/snic2020-6-155/INUPUT/DIR
OUT=/proj/snic2020-6-155/OUTPUT/DIR

# Create list of sample IDs by searching filenames
# Only need to run once!
find $DIR -maxdepth 1 -name '*.fastq.gz' -printf '%P\n' | awk -F '_' '{print $1}' | sort | uniq > $DIR/sample_IDs.txt

# Rename array IDs as sample IDs
SNAIL_ID=`sed -n "$SLURM_ARRAY_TASK_ID"p $DIR/sample_IDs.txt | awk '{print $1}'`

# Read trim call
fastp -i $DIR/${SNAIL_ID}_R1_001.fastq.gz \
        -I $DIR/${SNAIL_ID}_R2_001.fastq.gz \
        -o $OUT/${SNAIL_ID}_trimed_R1.fastq.gz \
        -O $OUT/${SNAIL_ID}_trimed_R2.fastq.gz \
        --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -h $OUT/${SNAIL_ID}_GeneWiz_fastp_summary_v2.html
### Note: fastp automatically identifies adapter sequences. This didn't work well with the trial data so I am specifying the sequences with --adapter_sequence flags.
