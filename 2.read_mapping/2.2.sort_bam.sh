#!/bin/bash

### Step 2.2 of short read processing pipeline
### Sorting reads and duplicate marking in BAM file. 
### This script was designed on a SBATCH system.

### James Reeve
### 2024-09-25

### Job parameters
#SBATCH -n 20
#SBATCH -t 2:00:00
#SBATCH -J Lsrp_tut_2.2
#SBATCH --output logfiles/Lsrp_sort_BAM.out
#SBATCH --error logfiles/Lsrp_sort_BAM.err

### Software
module load samtools/1.20

### File path
DIR=/path/to/your/bam

### Script
# '-@' assigns threads to samtools tasks
# i) Sort by name then position
samtools collate -@ 20 -Ou $DIR/CZA020_Ls_PE.bam |\
        samtools fixmate -@ 20 -m - - |\
        samtools sort -@ 20 - -o $DIR/CZA020_Ls_PE_sort.bam

# ii) Mark duplicates
samtools markdup -@ 20 -d 100 $DIR/CZA020_Ls_PE_sort.bam $DIR/CZA020_Ls_PE_DupMark.bam

# iii) Index and generate flagstat summary
samtools index -b $DIR/CZA020_Ls_PE_DupMark.bam
samtools flagstat $DIR/CZA020_Ls_PE_DupMark.bam > $DIR/CZA020_Ls_PE_DupMark.flagstat
