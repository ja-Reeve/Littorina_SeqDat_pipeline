#!/bin/bash

### This script maps trimmed WGS reads from Littorina saxatilis onto a reference
### genome freeze (De Jode, in prep). Multiple calls will be run in parallel using
### a sbatch array across all a list of sample names. 
### This script was designed to run on the Rackham server on UppMax.
### James Reeve - Göteborgs Universiter
### 29/03/2023

### Job parameters
#SBATCH -A snic2022-5-266
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 24:00:00
#SBATCH --array=1-137
#SBATCH -J read_map
#SBATCH --output=logfiles/read_map_%a.out

### Required software
module load bioinfo-tools
module load samtools/1.12
module load bwa/0.7.17

### The script:

# 1. Filepaths:
FASTQ=/proj/snic2020-6-155/TRIMMED/FASTQ/DIR
BAM=/proj/snic2020-6-155/OUTPUT/BAM/DIR
REF=/proj/snic2020-6-155/webexport/Lsax_new2023_ref

# 2. Assign sample names to array task IDs
SAMPLE=`sed -n "$SLURM_ARRAY_TASK_ID"p /SAMPLE/LIST/DIR/samples.txt | awk '{print $1}'`

# 3. Copy files to working drectory [optional]
# This avoids jumps in memory usage while navigating file paths
cp $REF/Lsax_genome_CLR_HiC_curated_freeze_1_2023_02_17.fasta $SNIC_TMP/
cp $FASTQ/${SNAIL_ID}*_R{1,2}.fastq.gz $SNIC_TMP/
cd $SNIC_TMP

# 4. Index the reference genome
bwa index Lsax_genome_CLR_HiC_curated_freeze_1_2023_02_17.fasta

# 5. Read mapping
# '-R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina"' not needed unless you plan to use GATK
bwa mem -M -t 20 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:illumina" \
        Lsax_genome_CLR_HiC_curated_freeze_1_2023_02_17.fasta \
        ${SAMPLE}_trimed_R1.fastq.gz ${SAMPLE}_trimed_R2.fastq.gz | \
        samtools view -b > ${SAMPLE}_PE.bam

# Error check [optional]
samtools view -H ${SAMPLE}_PE.bam | head -n 20

# 6. Sort and index the BAM
samtools sort -@ 20 -m 6400M -o $BAM/${SAMPLE}_PE_sorted.bam ${SAMPLE}_PE.bam
samtools index -b $BAM/${SAMPLE}_PE_sorted.bam
samtools flagstat $BAM/${SAMPLE}_PE_sorted.bam > $BAM/${SAMPLE}_PE.flagstat
