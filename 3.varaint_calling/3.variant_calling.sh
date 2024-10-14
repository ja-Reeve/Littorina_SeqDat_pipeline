#!/bin/bash

### Step 3 of short read processing pipeline
### Variant calling from a single BAM file. For multiple
### BAMs replace line 46 with "-b" and a list of BAM files.
### This script will take ages to run, 24h may not be long
### enough. Try using the code in "3.variant_calling_multithread.sh"
### to speed up your analysis.
### This script was designed on a SBATCH system.

### James Reeve
### 2024-10-02

### Job parameters
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -J Lsrp_tut_3
#SBATCH --output logfiles/Lsrp_variant_call.out
#SBATCH --error logfiles/Lsrp_variant_call.err

### Software
module load bcftools/1.20

### File path
DIR=/path/to/your/data
REF=/path/to/reference/genome

### Script
## i) Set parameters
# Chromosome [LG15]
CHR="CM074573.1"

# Calling
bcftools mpileup -Ou -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD \
        -f '"$REF"'/GCA_037325665.1_US_GU_Lsax_2.0_genomic.fna \
        -r $CHR \
        '"$DIR"'/CZA020_Ls_PE_DupMark.bam |\
        bcftools call -mOz -f GQ,GP -o '"$DIR"'/CZA020_Ls_$CHR.vcf.gz'
