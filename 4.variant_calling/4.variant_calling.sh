#!/bin/bash

### Step 4 of short read processing pipeline
### Filtering variants from a VCF. The thresholds used
### below are examples from the data used for this
### this tutorial. Follow step 4.4 "Choosing filtering 
### thresholds when filtering your own data. 
### This script was designed on a SBATCH system.

### James Reeve
### 2024-10-03

### Job parameters
#SBATCH -n 1
#SBATCH -t 20:00
#SBATCH -J Lsrp_tut_4
#SBATCH --output logfiles/Lsrp_variant_filter.out
#SBATCH --error logfiles/Lsrp_variant_filter.err

### Software
module load bcftools/1.20

### File path
DIR=/path/to/your/data

### Script
## i) Remove multiallelic sites, indels, and SNPs within 5bp
bcftools filter -Ou -g 5:indel,other $DIR/CZA020_Ls_CM074573.1.raw.vcf.gz |\
        bcftools view -Oz -M 2 -v snps > $DIR/CZA020_Ls_CM074573.1.noIndel.vcf.gz
# Add -m 2 to the previous line to remove invariant sites

# Count SNPs in filtered VCF
Nsnp=$(bcftools view -H $DIR/CZA020_Ls_CM074573.1.noIndel.vcf.gz | wc -l)
echo -e $Nsnp" SNPs after removing multiallelic sites, indels and nearby SNPs."


## ii) Remove sites with extremly high depth
bcftools filter -Oz -e 'INFO/DP > 40' $DIR/CZA020_Ls_CM074573.1.noIndel.vcf.gz > $DIR/CZA020_Ls_CM074573.1.depth.vcf.gz

# SNP count
Ndepth=$(bcftools view -H $DIR/CZA020_Ls_CM074573.1.depth.vcf.gz | wc -l)
echo -e $Ndpeth" SNPs after removing those with high total site depth [INFO/DP > 40]."


## iii) Remove sites with low mapping quality
bcftools filter -Oz -e 'MQ<30' $DIR/CZA020_Ls_CM074573.1.depth.vcf.gz > $DIR/CZA020_Ls_CM074573.1.MQ.vcf.gz

# SNP count
NMQ=$(bcftools view -H $DIR/CZA020_Ls_CM074573.1.MQ.vcf.gz | wc -l)
echo -e $NMQ" SNPs after removing those with low mapping quality [MQ < 30]."


## iv) Remove sites with low call quality
bcftools filter -Oz -e 'QUAL<30' $DIR/CZA020_Ls_CM074573.1.MQ.vcf.gz > $DIR/CZA020_Ls_CM074573.1.QUAL.vcf.gz

# SNP count
NQUAL=$(bcftools view -H $DIR/CZA020_Ls_CM074573.1.QUAL.vcf.gz | wc -l)
echo -e $NQUAL" SNPs after removing those with low call quality [QUAL < 30]."


## v) Remove sites with high strand bias
bcftools filter -Oz -e 'SP>3' $DIR/CZA020_Ls_CM074573.1.QUAL.vcf.gz > $DIR/CZA020_Ls_CM074573.1.SP.vcf.gz

# SNP count
NSP=$(bcftools view -H $DIR/CZA020_Ls_CM074573.1.SP.vcf.gz | wc -l)
echo -e $NSP" SNPs after removing those with high strand bias [SP > 3]." 


## vi) Set soft filter on genotypes
# Instead of removing sites, we will set genotypes
# which don't meet the threshold to missing (./.).
bcftools filter -Ou -S . -e 'FMT/GQ<20' $DIR/CZA020_Ls_CM074573.1.SP.vcf.gz |\
bcftools filter -Oz -S . -e 'FMT/DP<5' - > $DIR/CZA020_Ls_CM074573.1.soft_filt.vcf.gz
# No SNP count here, as this would not change from the last step.


## vii) Remove sites with more the 10% missing data
# Note: with only 1 sample the value is not important,
# but as the number of samples in the VCF goes up this
# threshold becomes relavent.
bcftools filter -Oz -e 'F_MISSING>0.1' $DIR/CZA020_Ls_CM074573.1.soft_filt.vcf.gz > $DIR/CZA020_Ls_CM074573.1.filtered.vcf.gz

# SNP count
Nfilt=$(bcftools view -H $DIR/CZA020_Ls_CM074573.1.filtered.vcf.gz | wc -l)
echo -e $Nfilt" SNPs after removing sites missing data [MISS > 10%]."
