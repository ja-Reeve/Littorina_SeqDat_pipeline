# Littorina short-read data processing pipeline
**James Reeve**

**2024-09-23**

A protocol for a standard pipeline for processing short-read seqeuncing data from *Littorina* snails. Each step was carefully considered to get the best results from short reads of a non-model species. The detailed desceription of this pipeline is published [here](https://www.protocols.io/view/a-standard-pipeline-for-processing-short-read-sequ-c6ygzftw).

This repository is intended as a tutorial for the protocol. All steps are tested on a single snail (CZA020_Ls) which was collected near the Tjärnö Marine Laboratory, Sweden (58.83091 N, 11.13305 E). Raw reads for this sample and the reference genome should be downloaded from NCBI before running the tutorial:

### Data access:
  Raw reads: NCBI-SRA [SAMN37254815](https://www.ncbi.nlm.nih.gov/biosample)
  
  Reference genome: NCBI-Genome [GCA_037325665.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_037325665.1/) 


### Protocol outline
![SNP_call_pipeline_Fig1](https://github.com/ja-Reeve/Littorina_SeqDat_pipeline/assets/82411887/a4f2633d-d8a3-47db-b244-2784c042aefc)
