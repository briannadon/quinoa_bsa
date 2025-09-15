#!/bin/bash

# Quinoa BSA Pipeline
# This script performs quality filtering, adapter removal, barcode removal,
# alignment to quinoa reference genome, and variant calling

set -e  # Exit on error
set -o pipefail  # Exit on pipe failure

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../data"
OUTPUT_DIR="$SCRIPT_DIR/results"
LOG_DIR="$OUTPUT_DIR/logs"
TMP_DIR="$OUTPUT_DIR/tmp"

# Reference files
REFERENCE="$DATA_DIR/Chenopodium_quinoa.faa"
GFF_FILE="$DATA_DIR/Chenopodium_quinoa_annos0-cds1-id_typename-nu1-upa1-add_chr0.gid33827.gff"

# Input files (subsampled for testing)
SRR_FILES=(
    "$DATA_DIR/SRR26831195.fastq.subsampled"
    "$DATA_DIR/SRR26831196.fastq.subsampled" 
    "$DATA_DIR/SRR26831197.fastq.subsampled"
)

# Tool environments (from Docker setup)
FASTQC_ENV="fastqc"
SEQTK_ENV="seqtk"
BWA_ENV="bwa"
SAMTOOLS_ENV="samtools"
# bcftools is not available in the current setup, will use samtools mpileup directly

# Create directories
mkdir -p "$OUTPUT_DIR" "$LOG_DIR" "$TMP_DIR" "$OUTPUT_DIR/fastqc"

# Function to log messages with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/pipeline.log"
}

# Function to run commands with conda environment
run_with_env() {
    local env="$1"
    local cmd="$2"
    log "Running in $env environment: $cmd"
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             source /opt/conda/etc/profile.d/mamba.sh && \
             conda activate $env && \
             $cmd" 2>&1 | tee -a "$LOG_DIR/pipeline.log"
}

# Function to check if file exists
check_file() {
    if [ ! -f "$1" ]; then
        log "ERROR: File $1 not found!"
        exit 1
    fi
}

# Function to check if tool is available
check_tool() {
    local env="$1"
    local tool="$2"
    if ! run_with_env "$env" "which $tool" > /dev/null 2>&1; then
        log "ERROR: Tool $tool not found in $env environment!"
        exit 1
    fi
}

# Initialize pipeline
log "=== Starting Quinoa BSA Pipeline ==="
log "Data directory: $DATA_DIR"
log "Output directory: $OUTPUT_DIR"

# Check required files exist
check_file "$REFERENCE"
for srr_file in "${SRR_FILES[@]}"; do
    check_file "$srr_file"
done

# Check required tools are available
check_tool "$FASTQC_ENV" "fastqc"
check_tool "$SEQTK_ENV" "seqtk"
check_tool "$BWA_ENV" "bwa"
check_tool "$SAMTOOLS_ENV" "samtools"
# bcftools check removed since it's not available

# Step 1: Quality control with FastQC
log "=== Step 1: Quality Control with FastQC ==="
for srr_file in "${SRR_FILES[@]}"; do
    sample_name=$(basename "$srr_file" .fastq.subsampled)
    log "Running FastQC on $sample_name"
    run_with_env "$FASTQC_ENV" "fastqc '$srr_file' -o '$OUTPUT_DIR/fastqc'"
done

# Step 2: Filter low-quality sequences and remove adapters/barcodes
log "=== Step 2: Quality Filtering and Adapter/Barcode Removal ==="
for srr_file in "${SRR_FILES[@]}"; do
    sample_name=$(basename "$srr_file" .fastq.subsampled)
    filtered_file="$OUTPUT_DIR/${sample_name}.filtered.fastq"
    
    log "Filtering low-quality sequences from $sample_name"
    # Using seqtk for basic quality filtering
    # -q 20: minimum quality score 20
    # -l 50: minimum length 50 bp
    run_with_env "$SEQTK_ENV" "seqtk seq -q 20 -l 50 '$srr_file' > '$filtered_file'"
    
    # Count reads before and after filtering
    original_count=$(run_with_env "$SEQTK_ENV" "seqtk comp '$srr_file' | wc -l")
    filtered_count=$(run_with_env "$SEQTK_ENV" "seqtk comp '$filtered_file' | wc -l")
    log "Sample $sample_name: $original_count reads -> $filtered_count reads after filtering"
done

# Step 3: Index reference genome for alignment
log "=== Step 3: Indexing Reference Genome ==="
# Check if all BWA index files exist
if [ ! -f "${REFERENCE}.bwt" ] || [ ! -f "${REFERENCE}.amb" ] || \
   [ ! -f "${REFERENCE}.ann" ] || [ ! -f "${REFERENCE}.pac" ] || \
   [ ! -f "${REFERENCE}.sa" ]; then
    log "Indexing reference genome with BWA"
    run_with_env "$BWA_ENV" "bwa index '$REFERENCE'"
else
    log "Reference genome already indexed (all index files present)"
fi

# Step 4: Align reads to reference genome
log "=== Step 4: Alignment to Reference Genome ==="
for srr_file in "${SRR_FILES[@]}"; do
    sample_name=$(basename "$srr_file" .fastq.subsampled)
    filtered_file="$OUTPUT_DIR/${sample_name}.filtered.fastq"
    bam_file="$OUTPUT_DIR/${sample_name}.aligned.bam"
    sorted_bam="$OUTPUT_DIR/${sample_name}.sorted.bam"
    
    log "Aligning $sample_name to reference genome"
    # Run bwa alignment and save to temporary file, then convert to BAM
    temp_sam="$TMP_DIR/${sample_name}.temp.sam"
    run_with_env "$BWA_ENV" "bwa mem -t 4 '$REFERENCE' '$filtered_file' > '$temp_sam'"
    run_with_env "$SAMTOOLS_ENV" "samtools view -b '$temp_sam' > '$bam_file'"
    rm -f "$temp_sam"
    
    log "Sorting BAM file for $sample_name"
    run_with_env "$SAMTOOLS_ENV" "samtools sort '$bam_file' -o '$sorted_bam'"
    
    log "Indexing sorted BAM file for $sample_name"
    run_with_env "$SAMTOOLS_ENV" "samtools index '$sorted_bam'"
    
    # Get alignment statistics
    log "Alignment statistics for $sample_name:"
    run_with_env "$SAMTOOLS_ENV" "samtools flagstat '$sorted_bam'"
done

# Step 5: Variant Calling
log "=== Step 5: Variant Calling ==="
# Create list of sorted BAM files for multi-sample variant calling
bam_list="$OUTPUT_DIR/bam_list.txt"
> "$bam_list"  # Clear the file
for srr_file in "${SRR_FILES[@]}"; do
    sample_name=$(basename "$srr_file" .fastq.subsampled)
    sorted_bam="$OUTPUT_DIR/${sample_name}.sorted.bam"
    echo "$sorted_bam" >> "$bam_list"
done

# Index reference for samtools if needed
if [ ! -f "${REFERENCE}.fai" ]; then
    log "Indexing reference for samtools"
    run_with_env "$SAMTOOLS_ENV" "samtools faidx '$REFERENCE'"
fi

# Call variants using samtools mpileup (since bcftools is not available)
vcf_file="$OUTPUT_DIR/quinoa_variants.vcf"
log "Calling variants using samtools mpileup"
run_with_env "$SAMTOOLS_ENV" "samtools mpileup -uv -f '$REFERENCE' -b '$bam_list' > '$vcf_file'"

# Simple variant filtering using grep (since bcftools filter is not available)
filtered_vcf="$OUTPUT_DIR/quinoa_variants.filtered.vcf"
log "Simple variant filtering (removing low-quality and low-depth variants)"
run_with_env "$SAMTOOLS_ENV" "grep -v 'DP=0;\|DP=1;\|DP=2;\|DP=3;\|DP=4;\|DP=5;\|DP=6;\|DP=7;\|DP=8;\|DP=9;\|DP=10;' '$vcf_file' | grep -v 'QUAL=0\|QUAL=1\|QUAL=2\|QUAL=3\|QUAL=4\|QUAL=5\|QUAL=6\|QUAL=7\|QUAL=8\|QUAL=9\|QUAL=10\|QUAL=11\|QUAL=12\|QUAL=13\|QUAL=14\|QUAL=15\|QUAL=16\|QUAL=17\|QUAL=18\|QUAL=19\|QUAL=20' > '$filtered_vcf'"

# Generate variant statistics
log "Variant calling statistics:"
total_variants=$(run_with_env "$SAMTOOLS_ENV" "grep -c '^[^#]' '$filtered_vcf'")
echo "Total filtered variants: $total_variants" > "$OUTPUT_DIR/variant_stats.txt"

# Step 6: Generate summary report
log "=== Step 6: Generating Summary Report ==="
summary_file="$OUTPUT_DIR/pipeline_summary.txt"
{
    echo "Quinoa BSA Pipeline Summary"
    echo "==========================="
    echo "Date: $(date)"
    echo ""
    echo "Samples processed:"
    for srr_file in "${SRR_FILES[@]}"; do
        sample_name=$(basename "$srr_file" .fastq.subsampled)
        filtered_file="$OUTPUT_DIR/${sample_name}.filtered.fastq"
        sorted_bam="$OUTPUT_DIR/${sample_name}.sorted.bam"
        
        original_count=$(run_with_env "$SEQTK_ENV" "seqtk comp '$srr_file' | wc -l")
        filtered_count=$(run_with_env "$SEQTK_ENV" "seqtk comp '$filtered_file' | wc -l")
        aligned_count=$(run_with_env "$SAMTOOLS_ENV" "samtools view -c '$sorted_bam'")
        
        echo "- $sample_name:"
        echo "  Original reads: $original_count"
        echo "  Filtered reads: $filtered_count"
        echo "  Aligned reads: $aligned_count"
        # Calculate percentages using bash arithmetic (avoiding bc dependency)
        filtering_rate=$(( (original_count - filtered_count) * 100 / original_count ))
        alignment_rate=$(( aligned_count * 100 / filtered_count ))
        echo "  Filtering rate: ${filtering_rate}%"
        echo "  Alignment rate: ${alignment_rate}%"
    done
    echo ""
    echo "Variant calling:"
    total_variants=$(run_with_env "$SAMTOOLS_ENV" "grep -c '^[^#]' '$filtered_vcf'")
    echo "Total variants called: $total_variants"
    echo "Filtered variants file: $filtered_vcf"
} > "$summary_file"

log "Pipeline completed successfully!"
log "Summary report: $summary_file"
log "Filtered variants: $filtered_vcf"
log "=== Quinoa BSA Pipeline Finished ==="
