#!/bin/bash

### Step 1.1 of Littorina short-read data processing pipeline:
### Read quality inspecition using fastqc and MultiQC.
### This script was designed on a SBATCH system.

### James Reeve
### 2024-09-23

### Job parameters
#SBATCH -n 4
#SBATCH -t 1:00:00
#SBATCH -J Lsrp_tut_1.1
#SBATCH --output logfiles/Lsrp_read_quality.out
#SBATCH --error logfiles/Lsrp_read_quality.err

### Software
module load fastqc/0.12.1
module load MultiQC/1.22.2

### File path
DIR=/path/to/testing/envrionment
FASTQ=/path/to/raw/fastq/data
TMP=/path/to/temporary/files

### Script
# fastqc
fastqc --dir $TMP -t 4 \
	$FASTQ/CZA020_Ls_*.fastq.gz \
	--outdir $DIR/fastqc

# MultiQC
multiqc $DIR/fastqc \
	-o $DIR/MultiQC \
	-n multiqc_report_CZA020_Ls.html
