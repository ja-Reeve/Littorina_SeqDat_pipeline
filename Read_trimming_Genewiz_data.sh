#!/bin/bash

### Trim whole genome sequence reads from Littorin saxatilis data.
### Reads were squenced by GeneWiz using an Illumina library
### preperation. Adaptor sequences were not clipped.
### James Reeve - Göteborgs Universitet
### 01/11/2022

### Job parameters
#SBATCH -A snic2022-5-266 
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 2:30:00
#SBATCH --array=1-137%20
#SBATCH -J BS_read_trim_GeneWiz
#SBATCH --output logfiles/BS_read_trim_GeneWiz_%a.log 

## Script
module load bioinfo-tools
module load trimmomatic/0.39
# File path
DIR=/proj/snic2020-6-155/James/fastq/Broad_scale/GeneWiz

# Create list of sample IDs by searching filenames
find $DIR -maxdepth 1 -name '*.fastq.gz' -printf '%P\n' | awk -F '_' '{print $1}' | sort | uniq > $DIR/sample_IDs.txt

# Rename array IDs as sample IDs
SNAIL_ID=`sed -n "$SLURM_ARRAY_TASK_ID"p $DIR/sample_IDs.txt | awk '{print $1}'`

# Read trim call
trimmomatic PE -threads 12 -phred33 \
        $DIR/${SNAIL_ID}_R1_001.fastq.gz $DIR/${SNAIL_ID}_R2_001.fastq.gz \
        -baseout $DIR/trimmed/${SNAIL_ID}.fastq.gz \
        ILLUMINACLIP:/proj/snic2020-6-155/James/ref_genome/Azenta_adapters.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:10:30 MINLEN:70
