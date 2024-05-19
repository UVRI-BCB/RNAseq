#!/bin/bash

# Define input files
reads_dir="/path/to/reads_directory"
index="/path/to/genome_index"
annotation="/path/to/annotation_file.gff3"

# Define output directory
output_dir="/path/to/output_directory"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Step 1: Quality control using FastQC
fastqc_dir="$output_dir/fastqc_output"
mkdir -p "$fastqc_dir"
for reads_file in "$reads_dir"/*_R1.fastq.gz; do
    reads_basename=$(basename "$reads_file" _R1.fastq.gz)
    fastqc "$reads_file" -o "$fastqc_dir"
done

# Step 2: Check if paired-end or single-end
paired_end=false
if ls "$reads_dir"/*_R2.fastq.gz 1> /dev/null 2>&1; then
    paired_end=true
    echo "Paired-end data detected."
else
    echo "Single-end data detected."
fi

# Step 3: Trimming using Trimmomatic
trimmed_dir="$output_dir/trimmed_reads"
mkdir -p "$trimmed_dir"
if [ "$paired_end" = true ]; then
    for reads_file in "$reads_dir"/*_R1.fastq.gz; do
        reads_basename=$(basename "$reads_file" _R1.fastq.gz)
        trimmomatic PE -phred33 "$reads_file" "${reads_file/_R1.fastq.gz/_R2.fastq.gz}" \
            "$trimmed_dir/${reads_basename}_R1_trimmed.fastq.gz" "$trimmed_dir/${reads_basename}_R1_unpaired.fastq.gz" \
            "$trimmed_dir/${reads_basename}_R2_trimmed.fastq.gz" "$trimmed_dir/${reads_basename}_R2_unpaired.fastq.gz" \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    done
else
    for reads_file in "$reads_dir"/*.fastq.gz; do
        reads_basename=$(basename "$reads_file" .fastq.gz)
        trimmomatic SE -phred33 "$reads_file" "$trimmed_dir/${reads_basename}_trimmed.fastq.gz" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    done
fi

# Step 4: Aligning to the human reference genome using STAR
aligned_dir="$output_dir/aligned_reads"
mkdir -p "$aligned_dir"
for trimmed_file in "$trimmed_dir"/*_R1_trimmed.fastq.gz; do
    trimmed_basename=$(basename "$trimmed_file" _R1_trimmed.fastq.gz)
    STAR --genomeDir "$index" --readFilesIn "$trimmed_file" "${trimmed_file/_R1_trimmed.fastq.gz/_R2_trimmed.fastq.gz}" --outFileNamePrefix "$aligned_dir/${trimmed_basename}_aligned_"
done

# Step 5: Generating gene counts using htseq-count
gene_counts_dir="$output_dir/gene_counts"
mkdir -p "$gene_counts_dir"
for aligned_file in "$aligned_dir"/*Aligned.sortedByCoord.out.bam; do
    aligned_basename=$(basename "$aligned_file" _Aligned.sortedByCoord.out.bam)
    htseq-count -f bam -s no -i gene_id "$aligned_file" "$annotation" > "$gene_counts_dir/${aligned_basename}_gene_counts.txt"
done
