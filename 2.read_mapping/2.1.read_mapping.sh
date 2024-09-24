#!/bin/bash

### Step 2.1 of short read processing pipeline
### Read mapping using BWA. This is a very time consuming step, 
### it may take over a day for the reads to map.
### This script was designed on a SBATCH system.

### James Reeve
### 2024-09-24

### Job parameters
#SBATCH -n 20
#SBATCH -t 2-0:00:00
#SBATCH -J Lsrp_tut_2.1
#SBATCH --output logfiles/Lsrp_read_map.out
#SBATCH --error logfiles/Lsrp_read_map.err

### Software
module load bwa/0.7.18
module load samtools/1.20

### File path
DIR=/path/to/testing/envrionment
REF=/path/to/reference/genome

### Script
# 1) Index reference genome
# Note: this only needs to be run once
bwa index $REF/GCA_037325665.1_US_GU_Lsax_2.0_genomic.fna.gz

# 2) Mapping
bwa mem \
        -M -t 20 \
        -R '@RG\tID:CZA020_Ls\tSM:CZA020_Ls\tPL:illumina' \
        $REF/GCA_037325665.1_US_GU_Lsax_2.0_genomic.fna.gz \
        $DIR/CZA020_Ls_R1.fastq.gz \
        $DIR/CZA020_Ls_R2.fastq.gz |\ # pipe to samtools to compress into BAM
samtools view -b > $DIR/CZA020_Ls_PE.bam
