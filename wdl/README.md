RNA-Seq Analysis Pipeline
This WDL (Workflow Description Language) script defines a pipeline for RNA-Seq analysis, including quality control, trimming, alignment, and gene expression quantification.

Requirements

* Cromwell
* FastQC
* MultiQC
* Trimmomatic
* STAR
* HTSeq

## Usage

* Provide your raw RNA-Seq reads in FASTQ format.
* Prepare the human reference genome index for alignment.
* Obtain the annotation file (in GFF3 format) for gene quantification.

Create a JSON file specifying the paths to your input files.

```bash
{
  "RNASeqPipeline.reads": "/path/to/reads",
  "RNASeqPipeline.index": "/path/to/genome_index",
  "RNASeqPipeline.annotation": "/path/to/annotation.gff3"
}
```

Run Cromwell: Execute the WDL script using Cromwell along with the input JSON file.

```
java -jar cromwell.jar run rnaseq_pipeline.wdl -i inputs.json
```
