version 1.0

# Define input file types
task FastQC {
  File input_fastq

  command {
    fastqc ${input_fastq}
  }

  output {
    File fastqc_output_html = "${input_fastq.basename}_fastqc.html"
    File fastqc_output_zip = "${input_fastq.basename}_fastqc.zip"
  }
}

task Trimmomatic {
  File input_fastq

  command {
    trimmomatic ...
  }

  output {
    File trimmed_fastq = "${input_fastq.basename}_trimmed.fastq.gz"
  }
}

task STARAlignment {
  File input_fastq
  File index

  command {
    STAR --genomeDir ${index} --readFilesIn ${input_fastq}
  }

  output {
    File alignment_bam = "${input_fastq.basename}.bam"
  }
}

task HTSeqCount {
  File input_bam
  File annotation

  command {
    htseq-count -f bam -s no -i gene_id ${input_bam} ${annotation} > gene_counts.txt
  }

  output {
    File gene_counts = "gene_counts.txt"
  }
}

# Define workflow
workflow RNASeqPipeline {
  File reads
  File index
  File annotation

  call FastQC { input: input_fastq = reads }
  call Trimmomatic { input: input_fastq = reads }
  call STARAlignment { input: input_fastq = Trimmomatic.trimmed_fastq, index = index }
  call HTSeqCount { input: input_bam = STARAlignment.alignment_bam, annotation = annotation }
}
