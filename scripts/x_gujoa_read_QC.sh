#!/bin/bash
echo "script start: download and initial sequencing read quality control"
date

# Navigate to the project directory
cd /proj/applied_bioinformatics/users/x_gujoa/MedBioinfo

# Create output directory if it doesn't exist
mkdir -p ./data/sra_fastq

# Download Illumina sequence data

# sqlite3 -batch -noheader -csv /proj/applied_bioinformatics/common_data/sample_collab.db "select run_accession from sample_annot spl left join sample2bioinformatician s2b using(patient_code) where username='x_gujoa';" > ./analyses/x_gujoa_run_accessions.txt

echo "Download fastq files"
srun --cpus-per-task=8 --time=00:30:00 \
	singularity exec /proj/applied_bioinformatics/common_data/meta.sif \
	xargs -a ./analyses/x_gujoa_run_accessions.txt -I {} \
	fastq-dump --split-files --gzip --readids --outdir ./data/sra_fastq/ --disable-multithreading {}

# manipulate raw sequencing FASTQ files with seqkit

echo "manipulate raw sequencing FASTQ files with seqkit"

# Step 1: Print statistics on each FASTQ file
echo "Printing statistics on each FASTQ file..."
srun --cpus-per-task=8 --time=00:30:00 singularity exec /proj/applied_bioinformatics/common_data/meta.sif \
	seqkit stats --threads 8 ./data/sra_fastq/*.fastq.gz > ./logs/seqkit_stats.log

# Step 2: Check for duplicate reads (dry-run mode)
echo "Checking for duplicate reads..."
srun --cpus-per-task=8 --time=00:30:00 singularity exec /proj/applied_bioinformatics/common_data/meta.sif \
	seqkit rmdup -n --threads 8 ./data/sra_fastq/*.fastq.gz > ./logs/seqkit_rmdup.log

# Step 3: Check for adapter sequences
# Adapter Read1: AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
# Adapter Read2: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

echo "Checking for full adapter sequences..."
srun --cpus-per-task=8 --time=00:30:00 singularity exec /proj/applied_bioinformatics/common_data/meta.sif \
	seqkit locate -i -p AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -p AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --threads 8 ./data/sra_fastq/*.fastq.gz > ./logs/seqkit_locate_full.log

echo "Checking for shortened adapter sequences..."
srun --cpus-per-task=8 --time=00:30:00 singularity exec /proj/applied_bioinformatics/common_data/meta.sif \
	seqkit locate -i -p AGATCGGAAGA -p GATCGGAAGAG --threads 8 ./data/sra_fastq/*.fastq.gz > ./logs/seqkit_locate_short.log

# Quality control the raw sequencing FASTQ files with fastQC
echo "Quality control the raw sequencing FASTQ files with fastQC"

mkdir -p ./analyses/fastqc

srun --cpus-per-task=2 --time=00:30:00 singularity exec /proj/applied_bioinformatics/common_data/meta.sif \
	xargs -I{} -a ./analyses/x_gujoa_run_accessions.txt fastqc -o ./analyses/fastqc --threads 2 data/sra_fastq/{}_1.fastq.gz data/sra_fastq/{}_2.fastq.gz


date
echo "script end."

