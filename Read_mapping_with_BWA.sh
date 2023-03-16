#!/bin/bash

### Second test run of read mapping of WGS onto new Littorina saxatilis genome.
### This will run the same command over multiple samples. Starting with just 
### Iberian samples, since there is a demand for these sequences. The memory
### usage is being optimised by copying files onto the working directory ($SNIC_TMP).
### James Reeve - Göteborgs Universitet
### 17/11/2022

### Job parameters
#SBATCH -A snic2022-5-266
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 18:00:00
#SBATCH --array=1-32
#SBATCH -J Test_read_map_on_Iberia
#SBATCH --output=logfiles/Test2_read_mapping_%a.out

### Required software
module load bioinfo-tools
module load samtools/1.12
module load bwa/0.7.17

### The script:  
# 1. Filepaths:
DIR=/proj/snic2020-6-155/James
FASTQ=/proj/snic2020-6-155/James/fastq/Broad_scale/GeneWiz/trimmed
BAM=/proj/snic2020-6-155/James/bam/Broad_scale/GeneWiz

# 2. Create list of Iberian samples:
# Note: 'SpE' = Spain East & 'PoS' = Portugal South
basename -a $FASTQ/{Sp,Po}*_1P.fastq.gz | awk -F'[_]' '{print $1}' > $DIR/Iberian_Snails.txt
SAMPLE=$DIR/Iberian_Snails.txt
SNAIL_ID=`sed -n "$SLURM_ARRAY_TASK_ID"p $SAMPLE | awk '{print $1}'`

# 3. Copy files to working drectory:
# This avoids jumps in memory usage while navigating file paths
cp $DIR/ref_genome/yahs.out_scaffolds_final_high_low_gc_contigs.fa* $SNIC_TMP/
cp $FASTQ/${SNAIL_ID}*.fastq.gz $SNIC_TMP/
cd $SNIC_TMP

# 3. Index reference genome:
# This test is run on a pre-release version on the new reference genome!
# Note: this line only needs to be run once!
#bwa index $DIR/ref_genome/yahs.out_scaffolds_final_high_low_gc_contigs.fa

# 4. Read alignment, then sorting and indexing each sample:
bwa mem -M -t 12 -R "@RG\tID:H757NDSX3CGAATACG+NTCGGTAA:4\tSM:UKNW-S5\tLB:Nextera\tPL:illumina" \
        yahs.out_scaffolds_final_high_low_gc_contigs.fa \
        ${SNAIL_ID}_1P.fastq.gz \
        ${SNAIL_ID}_2P.fastq.gz | \
        samtools view -b > ${SNAIL_ID}_PE.bam
        
# 4 Error check:
# See if all scaffolds are in BAM header
samtools view -H ${SNAIL_ID}_PE.bam | head -n 20

# 5. Wrangle the BAM file:
# Sort, index and summarise the BAM
# Also move output to storage project
samtools sort -@ 4 -m 10G -o $BAM/${SNAIL_ID}_PE_sorted.bam ${SNAIL_ID}_PE.bam
samtools index -b $BAM/${SNAIL_ID}_PE_sorted.bam
samtools flagstat $BAM/${SNAIL_ID}_PE_sorted.bam > $BAM/${SNAIL_ID}_PE_v2.flagstat
