version 1.0

# Define workflow inputs
workflow STAR_alignment {
    File reads
    File index
    File annotation
    String readType

    # Define workflow outputs
    output {
        Array[File] aligned_reads
        Array[File] gene_counts
    }

    # Define task for quality control using FastQC and MultiQC
    task qc {
        input {
            File reads
        }
        output {
            File html_report
        }
        command <<<
            fastqc -o fastqc_output ${reads}
            multiqc fastqc_output -o multiqc_output
        >>>
    }

    # Define task for trimming using Trimmomatic for paired-end data
    task trimmingPE {
        input {
            File reads
        }
        output {
            File[] trimmed_reads
            File[] unpaired_reads
        }
        command <<<
            trimmomatic PE -phred33 ${reads}/*_R1.fastq.gz ${reads}/*_R2.fastq.gz \
                trimmed_reads/${reads.baseName}_paired_R1.fastq.gz \
                trimmed_reads/${reads.baseName}_unpaired_R1.fastq.gz \
                trimmed_reads/${reads.baseName}_paired_R2.fastq.gz \
                trimmed_reads/${reads.baseName}_unpaired_R2.fastq.gz \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        >>>
    }

    # Define task for trimming using Trimmomatic for single-end data
    task trimmingSE {
        input {
            File reads
        }
        output {
            File trimmed_reads
        }
        command <<<
            trimmomatic SE -phred33 ${reads}/*_R1.fastq.gz \
                trimmed_reads/${reads.baseName}_trimmed.fastq.gz \
                LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        >>>
    }

    # Define task for aligning to the human reference genome using STAR for paired-end data
    task alignmentPE {
        input {
            File reads
            File index
        }
        output {
            File aligned_reads
            File gene_counts
        }
        command <<<
            STAR --genomeDir ${index} --readFilesIn ${reads}/*_paired_R1.fastq.gz ${reads}/*_paired_R2.fastq.gz --outFileNamePrefix aligned_reads/ --quantMode GeneCounts
        >>>
    }

    # Define task for aligning to the human reference genome using STAR for single-end data
    task alignmentSE {
        input {
            File reads
            File index
        }
        output {
            File aligned_reads
            File gene_counts
        }
        command <<<
            STAR --genomeDir ${index} --readFilesIn ${reads} --outFileNamePrefix aligned_reads/ --quantMode GeneCounts
        >>>
    }

    # Define workflow execution order
    scatter (read in [reads]) {
        if (readType == "PE") {
            call trimmingPE { input: reads = read }
            call alignmentPE { input: reads = trimmingPE.trimmed_reads, index = index }
        }
        else if (readType == "SE") {
            call trimmingSE { input: reads = read }
            call alignmentSE { input: reads = trimmingSE.trimmed_reads, index = index }
        }
    }

    aligned_reads = readType == "PE" ? alignmentPE.aligned_reads : alignmentSE.aligned_reads
    gene_counts = readType == "PE" ? alignmentPE.gene_counts : alignmentSE.gene_counts
}
