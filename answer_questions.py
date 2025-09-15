#!/usr/bin/env python3
"""
Script to answer questions from questions.yaml using results data
and update answers.yaml with the answers.
"""

import yaml
import re
import subprocess
import os
from pathlib import Path

def read_yaml(file_path):
    """Read YAML file and return data."""
    with open(file_path, 'r') as f:
        return yaml.safe_load(f)

def write_yaml(file_path, data):
    """Write data to YAML file."""
    with open(file_path, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)

def get_gc_content():
    """Get modal GC content from FastQC data files."""
    gc_values = []
    fastqc_dir = Path("submission/workflow/results/fastqc")
    
    for sample_dir in fastqc_dir.iterdir():
        if sample_dir.is_dir() and sample_dir.name.endswith("_fastqc"):
            data_file = sample_dir / "fastqc_data.txt"
            if data_file.exists():
                with open(data_file, 'r') as f:
                    content = f.read()
                    # Find %GC line
                    match = re.search(r'%GC\t(\d+)', content)
                    if match:
                        gc_values.append(int(match.group(1)))
    
    # Calculate average GC content across samples
    if gc_values:
        return round(sum(gc_values) / len(gc_values))
    return None

def get_mapping_stats():
    """Get mapping statistics from pipeline log."""
    log_file = Path("submission/workflow/results/logs/pipeline.log")
    mapping_stats = {}
    
    if log_file.exists():
        with open(log_file, 'r') as f:
            content = f.read()
            
            # Find mapping percentages for each sample
            samples = ["SRR26831195", "SRR26831196", "SRR26831197"]
            for sample in samples:
                # Find the mapped percentage line for each sample
                pattern = rf"{sample}.*mapped \((\d+\.\d+)% : N/A\)"
                match = re.search(pattern, content)
                if match:
                    mapping_stats[sample] = float(match.group(1))
    
    # Calculate average mapping percentage
    if mapping_stats:
        avg_mapped = sum(mapping_stats.values()) / len(mapping_stats)
        # Calculate unmapped percentage (100% - mapped %)
        avg_unmapped = 100 - avg_mapped
        return round(avg_mapped, 2), round(avg_unmapped, 2)
    
    return None, None

def count_chromosome4_variants():
    """Count variants on chromosome 4 (gi|4)."""
    vcf_file = Path("submission/workflow/results/quinoa_variants.filtered.vcf")
    count = 0
    
    if vcf_file.exists():
        # Use grep to count lines containing gi|4
        result = subprocess.run(
            ["grep", "-c", "^gi|4", str(vcf_file)],
            capture_output=True, text=True
        )
        if result.returncode == 0:
            count = int(result.stdout.strip())
    
    return count

def count_variants_with_alternate_alleles():
    """Count variant sites with at least 1 alternate allele call."""
    vcf_file = Path("submission/workflow/results/quinoa_variants.filtered.vcf")
    count = 0
    
    if vcf_file.exists():
        # Use grep to count lines that don't have <*> as the only alternate allele
        result = subprocess.run(
            ["grep", "-v", "<*>", str(vcf_file)],
            capture_output=True, text=True
        )
        if result.returncode == 0:
            # Count non-header lines that don't have <*>
            lines = result.stdout.strip().split('\n')
            count = len([line for line in lines if line and not line.startswith('#')])
    
    return count

def main():
    """Main function to answer questions and update answers.yaml."""
    print("Reading questions and current answers...")
    
    # Read questions and current answers
    questions = read_yaml("questions.yaml")
    answers = read_yaml("answers.yaml")
    
    # Get answers for each question
    print("Answering questions...")
    
    # Q1: Modal GC content
    if not answers['answers'].get('q1'):
        gc_content = get_gc_content()
        if gc_content is not None:
            answers['answers']['q1'] = gc_content
            print(f"Q1 answered: Modal GC content = {gc_content}%")
    
    # Q2: Percent of reads that don't map
    if not answers['answers'].get('q2'):
        mapped_percent, unmapped_percent = get_mapping_stats()
        if unmapped_percent is not None:
            answers['answers']['q2'] = round(unmapped_percent)
            print(f"Q2 answered: {unmapped_percent}% of reads don't map")
    
    # Q3: Percent of reads that uniquely map
    if not answers['answers'].get('q3'):
        mapped_percent, _ = get_mapping_stats()
        if mapped_percent is not None:
            answers['answers']['q3'] = mapped_percent
            print(f"Q3 answered: {mapped_percent}% of reads uniquely map")
    
    # Q4: Variants on chromosome 4
    if not answers['answers'].get('q4'):
        chr4_count = count_chromosome4_variants()
        if chr4_count > 0:
            answers['answers']['q4'] = chr4_count
            print(f"Q4 answered: {chr4_count} variants on chromosome 4")
    
    # Q5: Variant sites with at least 1 alternate allele
    if not answers['answers'].get('q5'):
        alt_allele_count = count_variants_with_alternate_alleles()
        if alt_allele_count > 0:
            answers['answers']['q5'] = alt_allele_count
            print(f"Q5 answered: {alt_allele_count} variant sites with alternate alleles")
    
    # Update answers.yaml
    print("Updating answers.yaml...")
    write_yaml("answers.yaml", answers)
    print("Answers updated successfully!")

if __name__ == "__main__":
    main()
