## RNA-Seq Analysis Pipelines

This repository contains pipelines for RNA-Seq analysis implemented in various languages and frameworks. Choose the pipeline implementation that best fits your requirements and infrastructure.

#### Required tools

- FastQC (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- MultiQC (https://multiqc.info)
- Trimmomatic (https://github.com/usadellab/Trimmomatic)
- STAR (https://github.com/alexdobin/STAR)
- HTSeq (https://htseq.readthedocs.io/)

#### The workflow 

![alt text](https://github.com/UVRI-BCB/RNAseq/blob/4f016e4e63c73002fb76cf3170a84a263305971e/images/rnaseq.png?raw=true)

#### Pipelines

- The [Nextflow pipeline](https://github.com/UVRI-BCB/RNAseq/tree/main/bash) is implemented using Nextflow, a domain-specific language for data-driven computational pipelines. It includes steps for quality control, trimming, alignment, and gene expression quantification.

- The [WDL pipeline](https://github.com/UVRI-BCB/RNAseq/tree/main/wdl) is implemented using Workflow Description Language (WDL), a language for describing data analysis workflows. It includes steps for quality control, trimming, alignment, and gene expression quantification.

- The [Bash pipeline](https://github.com/UVRI-BCB/RNAseq/tree/main/bash) is implemented as a bash script. It performs RNA-Seq analysis using native Linux tools and bash scripting. It includes steps for quality control, trimming, alignment, and gene expression quantification.



