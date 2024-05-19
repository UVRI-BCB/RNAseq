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

// Step 1: Check for Paired-end or Single-end data
process checkPairedEnd {
    input:
    file reads from reads_ch
    
    output:
    file 'paired_end_flag.txt' into flag_ch

    script:
    """
    if (new File(reads.getParent(), "*_R2.fastq.gz").list().size() > 0) {
        println "Paired-end data detected."
    } else {
        println "Single-end data detected."
    }
    """
}

// Step 2: Quality control using FastQC and MultiQC
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

// Step 3: Trimming using Trimmomatic for paired-end data
process trimmingPE {
    input:
    file(reads) from reads_ch

    output:
    file("trimmed_reads/*_paired_R1.fastq.gz") into trimmed_reads_ch
    file("trimmed_reads/*_unpaired_R1.fastq.gz") into unpaired_reads_ch

    script:
    """
    trimmomatic PE -phred33 $reads/*_R1.fastq.gz $reads/*_R2.fastq.gz \
        trimmed_reads/${reads.baseName}_paired_R1.fastq.gz \
        trimmed_reads/${reads.baseName}_unpaired_R1.fastq.gz \
        trimmed_reads/${reads.baseName}_paired_R2.fastq.gz \
        trimmed_reads/${reads.baseName}_unpaired_R2.fastq.gz \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Step 4: Trimming using Trimmomatic for single-end data
process trimmingSE {
    input:
    file(reads) from reads_ch

    output:
    file("trimmed_reads/*.fastq.gz") into trimmed_reads_ch

    script:
    """
    trimmomatic SE -phred33 $reads/*_R1.fastq.gz \
        trimmed_reads/${reads.baseName}_trimmed.fastq.gz \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// Step 5: Aligning to the human reference genome using STAR for paired-end data
process alignmentPE {
    input:
    file(reads) from trimmed_reads_ch
    file(index) from index_ch

    output:
    file("aligned_reads/*.bam") into alignment_ch

    script:
    """
    STAR --genomeDir $index --readFilesIn $reads/*_paired_R1.fastq.gz $reads/*_paired_R2.fastq.gz --outFileNamePrefix aligned_reads/ --quantMode GeneCounts
    """
}

// Step 6: Aligning to the human reference genome using STAR for single-end data
process alignmentSE {
    input:
    file(reads) from trimmed_reads_ch
    file(index) from index_ch

    output:
    file("aligned_reads/*.bam") into alignment_ch

    script:
    """
    STAR --genomeDir $index --readFilesIn $reads --outFileNamePrefix aligned_reads/ --quantMode GeneCounts
    """
}

// Define output channels
channel.fromPath('flag_ch')
channel.fromPath('fastqc_ch')
channel.fromPath('multiqc_ch')
channel.fromPath('trimmed_reads_ch')
channel.fromPath('alignment_ch')
