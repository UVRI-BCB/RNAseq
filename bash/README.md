## RNA-Seq Analysis Pipeline in Bash

This bash script performs RNA-Seq analysis, including quality control, trimming, alignment, and gene expression quantification.

### Requirements

* FastQC
* MultiQC
* Trimmomatic
* STAR
* HTSeq

### Usage

* Provide your raw RNA-Seq reads in FASTQ format.
* Prepare the human reference genome index for alignment.
* Obtain the annotation file (in GFF3 format) for gene quantification.
* Run the script

```
bash covac_rnaseq_pipeline.sh --reads /path/to/reads --index /path/to/index --annotation /path/to/annotation.gff3 --readType PE
```
