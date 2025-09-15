# Quinoa BSA Pipeline - Step-by-Step Documentation

## Overview
This pipeline processes SRR fastq reads through quality filtering, adapter/barcode removal, alignment to quinoa reference genome, and variant calling for bulk segregant analysis.

## Step 1: Quality Control with FastQC
**Purpose**: Assess initial read quality before processing
**Tool**: FastQC
**Action**: 
- Runs FastQC on all input SRR files
- Outputs HTML and ZIP reports to `results/fastqc/`
- Provides visual assessment of per-base quality, sequence duplication, adapter content, etc.

## Step 2: Quality Filtering and Adapter/Barcode Removal
**Purpose**: Remove low-quality sequences and technical sequences
**Tool**: seqtk
**Filtering Decisions**:
- **Minimum Quality Score**: 20 (Phred score ≥20, error rate ≤1%)
- **Minimum Read Length**: 50 bp (removes short fragments)
- **Adapter/Barcode Removal**: Basic filtering through quality trimming
**Action**:
- Filters each SRR file using `seqtk seq -q 20 -l 50`
- Counts reads before and after filtering to calculate filtering rate
- Outputs filtered FASTQ files to `results/`

## Step 3: Reference Genome Indexing
**Purpose**: Prepare reference genome for alignment
**Tool**: BWA
**Action**:
- Checks if BWA index files already exist (.amb, .ann, .bwt, .pac, .sa)
- **Smart skipping**: If all index files exist, skips indexing (saves ~39 minutes)
- If missing any index files, runs `bwa index` to create them
- Index files are stored alongside the reference FASTA

## Step 4: Alignment to Reference Genome
**Purpose**: Map filtered reads to quinoa reference genome
**Tools**: BWA MEM + samtools
**Alignment Parameters**:
- **Algorithm**: BWA MEM (for long reads)
- **Threads**: 4 (-t 4)
- **Output**: BAM format (binary, compressed)
**Processing Steps**:
1. BWA MEM alignment → temporary SAM file
2. samtools view converts SAM to BAM
3. samtools sort creates sorted BAM
4. samtools index creates BAM index
5. samtools flagstat generates alignment statistics
**Output**: Sorted, indexed BAM files for each sample

## Step 5: Variant Calling
**Purpose**: Identify genetic variants across samples
**Tool**: samtools mpileup (bcftools not available)
**Parameters**:
- `-uv`: Output uncompressed VCF, verbose
- `-f`: Reference genome
- `-b`: List of BAM files for multi-sample calling
**Variant Filtering Decisions**:
- **Minimum Depth**: 10 reads (DP > 10)
- **Minimum Quality**: 20 (QUAL > 20, Phred score ≥20)
- **Filtering Method**: grep-based filtering (since bcftools filter unavailable)
- Removes variants with DP=0-10 and QUAL=0-20

## Step 6: Summary Report Generation
**Purpose**: Provide comprehensive pipeline statistics
**Content**:
- Sample processing statistics (read counts)
- Filtering rates (percentage of reads removed)
- Alignment rates (percentage of filtered reads aligned)
- Variant calling results (total variants called)
- File locations for all outputs

## Key Filtering Decisions Summary
1. **Quality Threshold**: Phred score ≥20 (1% error rate)
2. **Length Threshold**: ≥50 bp minimum read length
3. **Variant Depth**: ≥10 supporting reads required
4. **Variant Quality**: Phred score ≥20 for called variants
5. **Multi-sample**: Variants called across all samples simultaneously

## Output Files
- `results/fastqc/`: Quality control reports
- `results/*.filtered.fastq`: Quality-filtered reads
- `results/*.sorted.bam`: Aligned, sorted BAM files
- `results/quinoa_variants.filtered.vcf`: Filtered variants
- `results/pipeline_summary.txt`: Comprehensive statistics
- `results/variant_stats.txt`: Variant calling metrics
- `results/logs/pipeline.log`: Complete execution log
