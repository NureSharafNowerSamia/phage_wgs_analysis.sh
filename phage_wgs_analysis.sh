#!/bin/bash

# Bacteriophage WGS analysis pipeline for Illumina NextSeq 550 paired-end reads

# Configurable variables
KRAKEN_DB="${KRAKEN_DB:-/path/to/default/db/minikrakendb}"
MASH_REF="REFSEQ_20240124_Bacteria_complete.msh"
THREADS=30 # Not suitable for smaller infrastructure
SPADES_MEMORY=120  # Not suitable for smaller infrastructure

# Function to check if a tool is installed
install_tool() {
    if ! command -v "$1" &> /dev/null
    then
        echo "$1 is not installed. Installing..."
        if [[ "$1" == "fastqc" ]]; then
            sudo apt update && sudo apt install -y fastqc
        elif [[ "$1" == "bbduk.sh" ]]; then
            sudo apt update && sudo apt install -y bbmap
        elif [[ "$1" == "fastp" ]]; then
            sudo apt update && sudo apt install -y fastp
        elif [[ "$1" == "spades.py" ]]; then
            sudo apt update && sudo apt install -y spades
        elif [[ "$1" == "reformat.sh" ]]; then
            sudo apt update && sudo apt install -y bbmap
        elif [[ "$1" == "quast.py" ]]; then
            sudo apt update && sudo apt install -y quast
        elif [[ "$1" == "bwa" ]]; then
            sudo apt update && sudo apt install -y bwa
        elif [[ "$1" == "samtools" ]]; then
            sudo apt update && sudo apt install -y samtools
        elif [[ "$1" == "kraken2" ]]; then
            sudo apt update && sudo apt install -y kraken2
        elif [[ "$1" == "bracken" ]]; then
            echo "Please download and install Bracken manually from https://ccb.jhu.edu/software/bracken/"
        elif [[ "$1" == "mash" ]]; then
            sudo apt update && sudo apt install -y mash
        else
            echo "No automated install available for $1. Please install it manually."
        fi
    else
        echo "$1 is already installed."
    fi
}

# List of required tools
tools=("fastqc" "bbduk.sh" "fastp" "spades.py" "reformat.sh" "quast.py" "bwa" "samtools" "kraken2" "bracken" "mash")

# Install each tool if necessary
for tool in "${tools[@]}"; do
    install_tool "$tool"
done

# Check Mash reference file
if [[ ! -f "$MASH_REF" ]]; then
    echo "Error: Mash reference file $MASH_REF not found. Please download it."
    exit 1
fi

# Check if all tools are installed
echo "All required tools installed. Running the main pipeline"

# Sample array
sam=("Sample_ID_1" "Sample_ID_2")

# Creating directories
mkdir -p report raw-data filtered-scaffolds mapping quast taxonomy mash

# Looping through each sample
for num in "${sam[@]}"
do
    # Check input files
    if [[ ! -f "${num}_R1.fastq.gz" || ! -f "${num}_R2.fastq.gz" ]]; then
        echo "Error: Input files for ${num} not found! Skipping..."
        continue
    fi

    # Running FastQC
    fastqc "${num}_R1.fastq.gz" "${num}_R2.fastq.gz" || { echo "FastQC failed for ${num}"; continue; }
    
    # Running BBDuk
    bbduk.sh in1="${num}_R1.fastq.gz" in2="${num}_R2.fastq.gz" \
             out1="${num}_cleaned_R1.fastq.gz" out2="${num}_cleaned_R2.fastq.gz" \
             threads="$THREADS" hdist=1 k=31 ref=phiX.fasta stats="${num}_stats.txt" || { echo "BBDuk failed for ${num}"; continue; }
    
    # Running Fastp
    fastp --in1 "${num}_cleaned_R1.fastq.gz" --in2 "${num}_cleaned_R2.fastq.gz" \
          --out1 "${num}_R1_trim.fastq.gz" --out2 "${num}_R2_trim.fastq.gz" \
          --json "${num}_fastp.json" --html "${num}_fastp.html" \
          --unpaired1 "${num}_R1_fail.fastq.gz" --unpaired2 "${num}_R2_fail.fastq.gz" \
          --thread 16 --detect_adapter_for_pe -r --cut_right_window_size 20 \
          --cut_right_mean_quality 30 -l 50 -g -5 20 -3 20 || { echo "Fastp failed for ${num}"; continue; }
    
    # Concatenate failed reads
    cat "${num}_R1_fail.fastq.gz" "${num}_R2_fail.fastq.gz" > "${num}_singles.fastq.gz"
    
    # Single report
    fastp --in1 "${num}_singles.fastq.gz" --out1 "${num}_single.fastq.gz" \
          --disable_adapter_trimming --json "${num}_single_fastp.json" \
          --html "${num}_single_fastp.html" --thread "$THREADS" || { echo "Fastp (singles) failed for ${num}"; continue; }
    
    # Remove intermediate files
    rm "${num}_cleaned_R1.fastq.gz" "${num}_cleaned_R2.fastq.gz" \
       "${num}_R1_fail.fastq.gz" "${num}_R2_fail.fastq.gz" "${num}_singles.fastq.gz"
    
    # Moving QC results to report
    mv "${num}_stats.txt" "${num}_fastp.json" "${num}_fastp.html" \
       "${num}_single_fastp.json" "${num}_single_fastp.html" \
       "${num}_R1_fastqc.html" "${num}_R1_fastqc.zip" "${num}_R2_fastqc.html" "${num}_R2_fastqc.zip" report/
    
    # Moving raw data files
    mv "${num}_R1.fastq.gz" "${num}_R2.fastq.gz" raw-data/
    
    # Assembly
    spades.py -t "$THREADS" -m "$SPADES_MEMORY" -s "${num}_single.fastq.gz" -1 "${num}_R1_trim.fastq.gz" -2 "${num}_R2_trim.fastq.gz" \
              --phred-offset 33 -o "assembly/${num}/" || { echo "SPAdes failed for ${num}"; continue; }
    
    # Reformat scaffolds
    reformat.sh in="assembly/${num}/scaffolds.fasta" out="filtered-scaffolds/${num}_scaffolds.fasta" minlength=500 || { echo "Reformat failed for ${num}"; continue; }
    
    # Quast analysis
    quast.py "filtered-scaffolds/${num}_scaffolds.fasta" -o "quast/${num}/" || { echo "QUAST failed for ${num}"; continue; }
    
    # Extract trimmed bases and single bases
    awk 'NR==18' "report/${num}_fastp.json" | tr -c -d 0-9 >> trimbases.txt
    awk 'NR==17' "report/${num}_single_fastp.json" | tr -c -d 0-9 >> singlebases.txt
    
    # Genome size extraction
    awk 'NR==10' "quast/${num}/report.txt" | tr -c -d 0-9 >> genomesize.txt
    
    # BWA indexing and mapping
    bwa index "filtered-scaffolds/${num}_scaffolds.fasta" || { echo "BWA index failed for ${num}"; continue; }
    bwa mem -t "$THREADS" "filtered-scaffolds/${num}_scaffolds.fasta" "${num}_R1_trim.fastq.gz" "${num}_R2_trim.fastq.gz" > "mapping/${num}.sam" || { echo "BWA mem failed for ${num}"; continue; }
    
    # Convert SAM to BAM, sort, and remove intermediates
    samtools view -b -S "mapping/${num}.sam" > "mapping/${num}.bam" || { echo "Samtools view failed for ${num}"; continue; }
    samtools sort "mapping/${num}.bam" -o "mapping/${num}_sorted.bam" || { echo "Samtools sort failed for ${num}"; continue; }
    rm "mapping/${num}.sam" "mapping/${num}.bam"
    
    # Taxonomy analysis with Kraken2 and Bracken
    kraken2 --use-names --threads "$THREADS" --db "$KRAKEN_DB" \
            --report "taxonomy/${num}_kraken.txt" \
            "filtered-scaffolds/${num}_scaffolds.fasta" > "taxonomy/${num}.kraken" || { echo "Kraken2 failed for ${num}"; continue; }
    
    bracken -d "$KRAKEN_DB" -i "taxonomy/${num}_kraken.txt" -l S \
            -o "taxonomy/${num}.bracken" || { echo "Bracken failed for ${num}"; continue; }
    
    # Mash analysis
    mash dist "filtered-scaffolds/${num}_scaffolds.fasta" "$MASH_REF" > "mash/${num}_mash.tab" || { echo "Mash failed for ${num}"; continue; }
done

echo "Pipeline completed!"
