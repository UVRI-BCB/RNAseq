nextflow.enable.dsl=2

params {
    // Define parameters
    file reads
    file index
    file annotation
    
    // Define channels
    Channel.fromPath(reads).set { reads_ch }
    Channel.fromPath(index).set { index_ch }
    Channel.fromPath(annotation).set { annotation_ch }
}

// Step 1: Quality control using FastQC and MultiQC
process qc {
    input:
    file(reads) from reads_ch

    output:
    file("fastqc_output/*.html") into fastqc_ch
    file("multiqc_output/multiqc_report.html") into multiqc_ch

    script:
    """
    fastqc -o fastqc_output $reads
    multiqc fastqc_output -o multiqc_output
    """
}

// Step 2: Trimming using Trimmomatic
process trimming {
    input:
    file(reads) from reads_ch

    output:
    file("trimmed_reads/*.fq.gz") into trimmed_reads_ch

    script:
    """
    trimmomatic ...
    """
}

// Step 3: Aligning to the human reference genome using STAR
process alignment {
    input:
    file(reads) from trimmed_reads_ch
    file(index) from index_ch

    output:
    file("aligned_reads/*.bam") into alignment_ch

    script:
    """
    STAR --genomeDir $index --readFilesIn $reads --outFileNamePrefix aligned_reads/
    """
}

// Step 4: Generating gene counts using htseq-count
process gene_counts {
    input:
    file(reads) from alignment_ch
    file(annotation) from annotation_ch

    output:
    file("gene_counts/*.txt") into gene_counts_ch

    script:
    """
    htseq-count -f bam -s no -i gene_id $reads $annotation > gene_counts/gene_counts.txt
    """
}

// Define output channels
channel.fromPath('fastqc_ch')
channel.fromPath('multiqc_ch')
channel.fromPath('trimmed_reads_ch')
channel.fromPath('alignment_ch')
channel.fromPath('gene_counts_ch')
