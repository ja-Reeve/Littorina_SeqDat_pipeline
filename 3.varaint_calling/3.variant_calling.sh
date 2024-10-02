#!/bin/bash

### Step 3 of short read processing pipeline
### Variant calling from a single BAM file. For multiple
### BAMs replace line 46 with "-b" and a list of BAM files.
### This script was designed on a SBATCH system.

### James Reeve
### 2024-10-02

### Job parameters
#SBATCH -n 20
#SBATCH -t 1:00:00
#SBATCH -J Lsrp_tut_3
#SBATCH --output logfiles/Lsrp_variant_call.out
#SBATCH --error logfiles/Lsrp_variant_call.err

### Software
module load bcftools/1.20
module load parallel/20230422

### File path
DIR=/path/to/your/data
REF=/path/to/reference/genome

### Script
## i) Set parameters
# Chromosome [LG15]
CHR="CM074573.1"
# Window size for parallel jobs
RegSize=1000000


## ii) create regions for each chromosome [this only needs to be run once]
cut -f1-2 $REF/GCA_037325665.1_US_GU_Lsax_2.0_genomic.fna.fai > $DIR/chrSize.txt
ChrL=$(grep -w "$CHR" $DIR/chrSize.txt | cut -f2)
paste -d '-' <(seq 1 $RegSize $ChrL) <(printf "$(seq $RegSize $RegSize $ChrL)\n$ChrL") > $DIR/${CHR}_regions.txt
sed -i -e 's/^/'"$CHR"':/' $DIR/${CHR}_regions.txt


## iii) Varinat call each window
cat $DIR/${CHR}_regions.txt | parallel -j20 '
        bcftools mpileup -Ou -a FORMAT/AD,FORMAT/DP,FORMAT/SP,INFO/AD \
                -f '"$REF"'/GCA_037325665.1_US_GU_Lsax_2.0_genomic.fna \
                -r ${} \
                '"$DIR"'/CZA020_Ls_PE_DupMark.bam |\
        bcftools call -mOz -f GQ,GP -o '"$DIR"'/CZA020_Ls_{}.vcf.gz'
# Note: output VCFs will be named by region (e.g., CZA020_Ls_chr1:1-1000000.vcf.gz


## iv) Concatenate VCFs of each region
ls $DIR/CZA020_Ls_$CHR:*.vcf.gz > $DIR/list.of.vcf.windows.txt
bcftools concat -f $DIR/list.of.vcf.windows.txt |\
bcftools sort -Oz - > $DIR/CZA020_Ls_$CHR.raw.vcf.gz


## Clean up intermediate files
rm $DIR/CZA020_Ls_$CHR:*.vcf.gz $DIR/list.of.vcf.windows.txt $DIR/${CHR}_regions.txt
