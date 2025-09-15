#!/bin/bash

# Script to generate summary report from completed pipeline results
# This can be run independently after the main pipeline completes

set -e  # Exit on error

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/../data"
OUTPUT_DIR="$SCRIPT_DIR/results"
LOG_DIR="$OUTPUT_DIR/logs"

# Reference files
REFERENCE="$DATA_DIR/Chenopodium_quinoa.faa"

# Input files
SRR_FILES=(
    "$DATA_DIR/SRR26831195.fastq.subsampled"
    "$DATA_DIR/SRR26831196.fastq.subsampled" 
    "$DATA_DIR/SRR26831197.fastq.subsampled"
)

# Tool environments
SEQTK_ENV="seqtk"
SAMTOOLS_ENV="samtools"

# Function to log messages with timestamp
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/summary.log"
}

# Function to run commands with conda environment
run_with_env() {
    local env="$1"
    local cmd="$2"
    log "Running in $env environment: $cmd"
    bash -c "source /opt/conda/etc/profile.d/conda.sh && \
             source /opt/conda/etc/profile.d/mamba.sh && \
             conda activate $env && \
             $cmd" 2>&1 | tee -a "$LOG_DIR/summary.log"
}

# Initialize summary generation
log "=== Generating Summary Report from Completed Results ==="
log "Output directory: $OUTPUT_DIR"

# Check if required files exist
if [ ! -d "$OUTPUT_DIR" ]; then
    log "ERROR: Output directory $OUTPUT_DIR not found!"
    exit 1
fi

# Step 6: Generate summary report
log "=== Generating Summary Report ==="
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
        
        # Get counts from existing files
        original_count=$(run_with_env "$SEQTK_ENV" "seqtk comp '$srr_file' | wc -l")
        filtered_count=$(run_with_env "$SEQTK_ENV" "seqtk comp '$filtered_file' | wc -l")
        aligned_count=$(run_with_env "$SAMTOOLS_ENV" "samtools view -c '$sorted_bam' 2>/dev/null || echo 0")
        
        echo "- $sample_name:"
        echo "  Original reads: $original_count"
        echo "  Filtered reads: $filtered_count"
        echo "  Aligned reads: $aligned_count"
        # Calculate percentages using bash arithmetic
        filtering_rate=$(( (original_count - filtered_count) * 100 / original_count ))
        alignment_rate=$(( aligned_count * 100 / filtered_count ))
        echo "  Filtering rate: ${filtering_rate}%"
        echo "  Alignment rate: ${alignment_rate}%"
    done
    echo ""
    echo "Variant calling:"
    filtered_vcf="$OUTPUT_DIR/quinoa_variants.filtered.vcf"
    if [ -f "$filtered_vcf" ]; then
        total_variants=$(run_with_env "$SAMTOOLS_ENV" "grep -c '^[^#]' '$filtered_vcf'")
        echo "Total variants called: $total_variants"
        echo "Filtered variants file: $filtered_vcf"
    else
        echo "Variant calling not completed yet"
    fi
} > "$summary_file"

log "Summary report generated successfully!"
log "Summary file: $summary_file"
log "=== Summary Generation Finished ==="
