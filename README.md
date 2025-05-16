# Bacteriophage WGS Analysis Pipeline

## Description
This Bash pipeline processes paired-end Illumina NextSeq 550 FASTQ files (R1 and R2) for bacteriophage whole-genome sequencing (WGS) analysis. It performs quality control, adapter trimming, assembly, quality assessment, read mapping, taxonomic classification, and genomic distance analysis.

## Prerequisites
- **Operating System**: Linux (Ubuntu recommended)
- **Dependencies**:
  - FastQC
  - BBMap (includes BBDuk, reformat.sh)
  - Fastp
  - SPAdes
  - QUAST
  - BWA
  - Samtools
  - Kraken2
  - Bracken (manual installation may be required)
  - Mash
- **Databases**:
  - Kraken2 database (e.g., MiniKraken)
  - Mash reference (`REFSEQ_20240124_Bacteria_complete.msh`)
  - phiX.fasta (for adapter trimming)


## Input
- Paired-end FASTQ files (e.g., `sample_R1.fastq.gz`, `sample_R2.fastq.gz`) in the working directory.
- Edit the `sam` array in the script to list your sample prefixes (e.g., `sam=("sample1" "sample2")`).

## Usage
1. Place FASTQ files in the working directory.
2. Update the `sam` array in `phage_wgs_pipeline.sh` with your sample names.
3. Run the pipeline:
   ```bash
   bash phage_wgs_pipeline.sh
   ```
4. Outputs are organized in directories:
   - `report/`: QC reports (FastQC, Fastp, BBDuk stats)
   - `raw-data/`: Original FASTQ files
   - `filtered-scaffolds/`: Assembled scaffolds (>500 bp)
   - `mapping/`: Sorted BAM files
   - `quast/`: Assembly quality metrics
   - `taxonomy/`: Kraken2 and Bracken taxonomy reports
   - `mash/`: Mash distance results
