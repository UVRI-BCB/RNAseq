# RNAseq

RNA-Seq Analysis Pipeline in Nextflow

This Nextflow script executes a complete RNA-Seq analysis pipeline, including quality control, trimming, alignment, and gene expression quantification.

Requirements:
- Nextflow (https://www.nextflow.io/)
- Trimmomatic (https://github.com/usadellab/Trimmomatic)
- STAR (https://github.com/alexdobin/STAR)
- HTSeq (https://htseq.readthedocs.io/en/release_0.11.1/)

Usage:
1. Prepare your input files:
   - Provide your raw RNA-Seq reads in FASTQ format.
   - Prepare the human reference genome index for alignment.
   - Obtain the annotation file (in GFF3 format) for gene quantification.

2. Edit the parameters in the `rnaseq_pipeline.nf` script:
   - Set the paths to your input files in the `params` block.
   - Update any command options if necessary.

3. Run the pipeline:

```bash
nextflow run covac_rnaseq_pipeline.nf --reads "/path/to/reads.fastq.gz" --index "/path/to/genome_index" --annotation "/path/to/annotation.gff3"
```
