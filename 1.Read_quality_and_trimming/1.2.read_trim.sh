#!/bin/bash

### Step 1.2 of short read processing pipeline
### Read trimming using fastp.
### This script was designed on a SBATCH system.

### James Reeve
### 2024-09-23

### Job parameters
#SBATCH -n 4
#SBATCH -t 1:00:00
#SBATCH -J Lsrp_tut_1.2
#SBATCH --output logfiles/Lsrp_read_trim.out
#SBATCH --error logfiles/Lsrp_read_trim.err

### Software
module load UPPMAX/1.0.0 bioinfo-tools
module load fastp/0.23.4

### File path
DIR=/path/to/testing/envrionment
FASTQ=/path/to/raw/fastq/data

### Script
fastp \
        -i $FASTQ/CZA020_Ls_merged_1.fastq.gz \     # Input forward read
        -I $FASTQ/CZA020_Ls_merged_2.fastq.gz \     # Input reverse read
        -o $DIR/CZA020_Ls_R1.fastq.gz \             # Output forward read
        -O $DIR/CZA020_Ls_R2.fastq.gz \             # Output reverse read
        --thread 4 -g -c -y 30 \
        --html $DIR/trim_quality_report.html \
        --json $DIR/trim_quality_report.json \
        --report_title CZA020_Ls_trim_report
