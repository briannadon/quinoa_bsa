

# Quinoa Bulk Segregant Analysis

## 1. Data sources

This task is based on publicly available sequencing data from a study of **Quinoa bulks segregant analysis**. The dataset includes RNAseq samples and was originally sequenced using **Illumina sequencing platform**.

The FASTQ files are stored in `data/` and are used as the inputs for the workflow. They are subsampled down to 1M reads via `seqtk sample -s100 sample.fq.gz 1000000`.

---

## 2. How to download

The data was downloaded from NCBI SRA using the SRA Toolkit. The SRA accession numbers are listed in `SraAccList.csv`.

### Example using SRA Toolkit

```bash
# Install SRA Toolkit
conda install -c bioconda sra-tools

# Download specific SRA accession
prefetch SRR24303833
fasterq-dump SRR24303833 --split-files
```

## 3. Running the workflow

The workflow is executed using submission/workflow/pipeline.sh. There is also a script that will run it in a pre-built docker container if needed (see Dockerfile for more)

