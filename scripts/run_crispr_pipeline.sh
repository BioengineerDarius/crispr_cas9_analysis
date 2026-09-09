#!/usr/bin/env bash
set -euo pipefail

echo "=================================================="
echo " Starting Automated CRISPR-Cas9 Target Pipeline  "
echo "=================================================="

# Ensure output directories exist
mkdir -p results/sgRNA_design results/off_target_analysis qc data

# 1. Target Sequence PAM Scanning & Tm Calculation
echo "[Step 1/4] Scanning target genome for 5'-NGG-3' PAM sites..."
grep -E -o "[ATCG]{20}[ATCG]GG" data/target_genome.fasta > results/sgRNA_design/candidate_sgRNAs.txt

awk '{
    g=gsub(/[GCgc]/,"",$1); 
    a=gsub(/[ATat]/,"",$1); 
    print $1, "Tm:" (4*g)+(2*a)"C"
}' results/sgRNA_design/candidate_sgRNAs.txt > results/sgRNA_design/candidate_tm_scores.txt

# 2. BLAST Database & Off-Target Profiling
echo "[Step 2/4] Executing BLAST off-target screening..."
awk '{print ">sgRNA_"NR"\n"$1}' results/sgRNA_design/candidate_sgRNAs.txt > results/sgRNA_design/candidates.fasta

if [ ! -f "data/ref_genome_db.nhr" ]; then
    makeblastdb -in data/target_genome.fasta -dbtype nucl -out data/ref_genome_db
fi

blastn -query results/sgRNA_design/candidates.fasta -db data/ref_genome_db -word_size 7 -evalue 1000 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -out results/off_target_analysis/off_target_matches.tsv

# 3. Read QC and Trimming
echo "[Step 3/4] Trimming sequencing reads with fastp..."
fastp -i data/SRR13834165_1.fastq -I data/SRR13834165_2.fastq -o data/trimmed_R1.fq.gz -O data/trimmed_R2.fq.gz --qualified_quality_phred 20 --html qc/fastp_report.html

# 4. Alignment & Indexing
echo "[Step 4/4] Aligning trimmed reads with BWA and indexing BAM..."
if [ ! -f "data/target_genome.fasta.bwt" ]; then
    bwa index data/target_genome.fasta
fi

bwa mem -t 4 data/target_genome.fasta data/trimmed_R1.fq.gz data/trimmed_R2.fq.gz | samtools view -bS - | samtools sort -o results/aligned_cleavage.bam
samtools index results/aligned_cleavage.bam

echo "=================================================="
echo "    CRISPR-Cas9 Pipeline Execution Complete!     "
echo "=================================================="

