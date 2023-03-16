#!/bin/bash

### Testing DeepVariant caller on Littorina saxatilis data
### Website: https://github.com/google/deepvariant

### Trial 1.2: single call on all individuals at one contig.
### This follows the same regions as Octopus Trial 1.8, so
### that the varaint callers can be directly compared.
### Edit [2022-12-12]: using bcftools to merge DeepVariant
### calls

### James Reeve - Göteborgs Universitet
### 07/12/2022

### Job parameters
#SBATCH -A snic2022-5-266
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 1-00:00:00
#SBATCH -J DeepVar_Trial_1.2
#SBATCH --output logfiles/DeepVariant_Trial1.2_v2.log 

### Required software
module load bioinfo-tools
module load DeepVariant/1.3.0
module load bcftools/1.14

### Command
DIR=/proj/snic2020-6-155/James
BAMLIST=$DIR/bam_list_SeanWGS.txt

### Run DeepVariant for each individual
while read BAM
do
        deepvariant --model_type=WGS \
                --ref=$DIR/ref_genome/Littorina_scaffolded_PacBio_run2_7_Oct_2016_unmasked.fasta \
                --reads=$DIR/bam/Seans_bams/${BAM}.RG.bam \
                --regions "Contig45579:0-453,697" \
                --output_vcf=$DIR/vcf/DeepVariant/DeepVariant_trial1.2_${BAM}.vcf \
                --output_gvcf=$DIR/vcf/DeepVariant/DeepVariant_trial1.2_${BAM}.vcf.gz \
                --num_shards=10

        rm $DIR/vcf/DeepVariant/DeepVariant_trial1.2_*.vcf
done < <(awk -F. '{print $1}' $BAMLIST)

### Make list of VCF files
ls $DIR/vcf/DeepVariant/DeepVariant_trial1*.vcf.gz > $DIR/vcf/DeepVariant/Trial1.2_vcf_list.txt

### Merge vcfs with bcftools
bcftools merge -m both --threads 10 -Oz \
        -o $DIR/vcf/DeepVariant/DeepVariant_trial1.2.vcf.gz \
        -l $DIR/vcf/DeepVariant/Trial1.2_vcf_list.txt
