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

# Merging paired end reads

echo "Merging paired end reads"

mkdir -p ./data/merged_pairs

srun --cpus-per-task=2 --time=01:00:00 singularity exec /proj/applied_bioinformatics/common_data/meta.sif \
	xargs -a ./analyses/x_gujoa_run_accessions.txt -I{} -n 1 bash -c 'flash -t 2 -z -o {}.flash -d ./data/merged_pairs ./data/sra_fastq/{}_1.fastq.gz ./data/sra_fastq/{}_2.fastq.gz 2>&1 | tee -a ./analyses/x_gujoa_flash.log'

# Read mapping
## Check for PhiX contamination (and more...)

echo "Check for PhiX contamination (and more...)"

mkdir -p ./data/reference_seqs

# download sequence data from the NCBI leveraging the ncbi edirect tool kit:
singularity exec /proj/applied_bioinformatics/common_data/meta.sif efetch -db nuccore -id NC_001422 -format fasta > ./data/reference_seqs/PhiX_NC_001422.fna

## Build the Bowtie2 index from the PhiX genome sequence:
### Note: such index creation step only ever needs to be run once for any given set of reference sequences
"echo Build the Bowtie2 index from the PhiX genome sequence:"

mkdir -p ./data/bowtie2_DBs

run singularity exec /proj/applied_bioinformatics/common_data/meta.sif bowtie2-build -f ./data/reference_seqs/PhiX_NC_001422.fna ./data/bowtie2_DBs/PhiX_bowtie2_DB

# Align Merged Reads to PhiX Genome
echo "Align Merged Reads to PhiX Genome"

srun --cpus-per-task=8 singularity exec /proj/applied_bioinformatics/common_data/meta.sif \
	bowtie2 -x ./data/bowtie2_DBs/PhiX_bowtie2_DB -U ./data/merged_pairs/*.extendedFrags.fastq.gz -S ./analyses/bowtie/x_gujoa_merged2PhiX.sam --threads 8 --no-unal 2>&1 | tee ./analyses/bowtie/x_gujoa_bowtie_merged2PhiX.log

# Align Reads to SARS-CoV-2 Genome
echo "Align Reads to SARS-CoV-2 Genome"
echo "Download the SARS-CoV-2 Genome sequence"

singularity exec /proj/applied_bioinformatics/common_data/meta.sif efetch -db nuccore -id NC_045512 -format fasta > ./data/reference_seqs/SC2_NC_045512.fna

echo "Build the Bowtie2 index for SARS-CoV-2:"

srun singularity exec /proj/applied_bioinformatics/common_data/meta.sif bowtie2-build -f ./data/reference_seqs/SC2_NC_045512.fna ./data/bowtie2_DBs/SC2_bowtie2_DB

echo "Align merged reads to the SARS-CoV-2 genome:"
srun --cpus-per-task=8 singularity exec /proj/applied_bioinformatics/common_data/meta.sif \
	bowtie2 -x ./data/bowtie2_DBs/SC2_bowtie2_DB -U ./data/merged_pairs/*.extendedFrags.fastq.gz -S ./analyses/bowtie/x_gujoa_merged2SC2.sam --threads 8 --no-unal 2>&1 | tee ./analyses/bowtie/x_gujoa_bowtie_merged2SC2.log


date
echo "script end."

