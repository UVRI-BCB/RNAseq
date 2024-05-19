#!/bin/bash

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --reads)
            reads_dir=$2
            shift 2
            ;;
        --index)
            index=$2
            shift 2
            ;;
        --annotation)
            annotation=$2
            shift 2
            ;;
        --readType)
            readType=$2
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Check if all required parameters are provided
if [[ -z $reads_dir || -z $index || -z $annotation || -z $readType ]]; then
    echo "Usage: $0 --reads <reads_directory> --index <genome_index> --annotation <annotation_file.gff3> --readType <SE or PE>"
    exit 1
fi

# Create output directories
mkdir -p fastqc_output multiqc_output trimmed_reads aligned_reads gene_counts

# Step 1: Quality control using FastQC and MultiQC
echo "Performing quality control..."
fastqc -o fastqc_output $reads_dir/*.fastq.gz
multiqc fastqc_output -o multiqc_output

# Step 2: Trimming using Trimmomatic
echo "Performing trimming..."
if [ "$readType" == "PE" ]; then
    for reads_file in $reads_dir/*_R1.fastq.gz; do
        reads_basename=$(basename "$reads_file" _R1.fastq.gz)
        trimmomatic PE -phred33 "$reads_file" "${reads_file/_R1.fastq.gz/_R2.fastq.gz}" \
            trimmed_reads/"${reads_basename}_paired_R1.fastq.gz" unpaired_reads/"${reads_basename}_unpaired_R1.fastq.gz" \
            trimmed_reads/"${reads_basename}_paired_R2.fastq.gz" unpaired_reads/"${reads_basename}_unpaired_R2.fastq.gz" \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    done
else
    for reads_file in $reads_dir/*.fastq.gz; do
        reads_basename=$(basename "$reads_file" .fastq.gz)
        trimmomatic SE -phred33 "$reads_file" \
            trimmed_reads/"${reads_basename}_trimmed.fastq.gz" \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    done
fi

# Step 3: Aligning to the human reference genome using STAR
echo "Performing alignment..."
if [ "$readType" == "PE" ]; then
    for trimmed_file in trimmed_reads/*_paired_R1.fastq.gz; do
        trimmed_basename=$(basename "$trimmed_file" _paired_R1.fastq.gz)
        STAR --genomeDir "$index" --readFilesIn "$trimmed_file" "${trimmed_file/_R1.fastq.gz/_R2.fastq.gz}" \
            --outFileNamePrefix aligned_reads/${trimmed_basename}_aligned_ --quantMode GeneCounts
    done
else
    for trimmed_file in trimmed_reads/*.fastq.gz; do
        trimmed_basename=$(basename "$trimmed_file" _trimmed.fastq.gz)
        STAR --genomeDir "$index" --readFilesIn "$trimmed_file" \
            --outFileNamePrefix aligned_reads/${trimmed_basename}_aligned_ --quantMode GeneCounts
    done
fi

# Step 4: Generating gene counts using htseq-count
echo "Generating gene counts..."
for aligned_file in aligned_reads/*ReadsPerGene.out.tab; do
    aligned_basename=$(basename "$aligned_file" _ReadsPerGene.out.tab)
    htseq-count -f tab -s no -i gene_id "$aligned_file" "$annotation" > gene_counts/${aligned_basename}_gene_counts.txt
done

echo "Workflow completed!"
