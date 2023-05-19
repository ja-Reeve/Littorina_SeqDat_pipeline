#!/bin/bash

### Script to identify and remove duplicated reads from WGS from Littorina saxatilis.
### Multiple samples are run in parallel on separate cores using a sbatch array. If you
### do not want to remove duplicates, take out the `-r`option in the markdups command.

### James Reeve - Göteborgs Universitet
### 19/05/2023

### Job parameters
#SBATCH -A snic2022-5-266
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 8:00:00
#SBATCH -J BS_markdup
#SBATCH --array=1-137
#SBATCH --output=logfiles/markdups_%a.out

### Required software
module load bioinfo-tools
module load samtools/1.12

### The script:

# 1. Filepaths:
BAM=/path/to/bams

# 2. Change array ID to sample ID
ID=`sed -n "$SLURM_ARRAY_TASK_ID"p /path/to/sample_IDs.txt | awk '{print $1}'`

# 2. fixmate and re-sort BAMs
# Sort reads by name
samtools sort -n -o $TMPDIR/${ID}_PE_sort1.bam $BAM/${ID}_PE_sorted.bam
# Add read pair details
samtools fixmate -m $TMPDIR/${ID}_PE_sort1.bam $TMPDIR/${ID}_PE_fixmated.bam
# Sort reads by position
samtools sort -o $TMPDIR/${ID}_PE_sort2.bam $TMPDIR/${ID}_PE_fixmated.bam

# 3. samtools-markdup
samtools markdup -d 100 $TMPDIR/${ID}_PE_sort2.bam $BAM/${ID}_PE_noDup.bam

# 4. Index and generate flagstats
samtools index -b $BAM/${ID}_PE_noDup.bam
samtools flagstat $BAM/${ID}_PE_noDup.bam > $BAM/${ID}_PE_noDup.flagstat
