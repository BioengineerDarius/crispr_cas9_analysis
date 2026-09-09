# CRISPR-Cas9 Target Identification & Off-Target Profiling Pipeline

An production-grade bioinformatics workflow for discovering single-guide RNA (sgRNA) targets, scanning for canonical $5'-NGG-3'$ PAM motifs, evaluating off-target alignment risk using
 BLAST, and processing high-throughput sequencing reads to quantify Cas9 cleavage.

## Workflow Overview
1. **Target Discovery & PAM Scanning**: Regular expression scanning for 20-nt protospacers preceding $5'-NGG-3'$ PAM motifs with thermodynamic $T_m$ calculation.
2. **Off-Target Assessment**: Local `blastn` short-read alignment against reference host genome databases.
3. **Read Quality Control**: `fastp` automated adapter trimming and quality filtering ($Q \ge 20$).
4. **Alignment & Indexing**: High-throughput read mapping via `BWA-MEM` and `SAMtools` BAM coordinate sorting/indexing.

## Repository Architecture
```text
crispr_cas9_analysis/
├── results/
│   ├── sgRNA_design/          # Candidate sgRNA sequences & Tm scores
│   └── off_target_analysis/   # BLAST hit alignment tables
├── scripts/
│   └── run_crispr_pipeline.sh # Master pipeline driver script
├── .gitignore                 # Version control exclusion rules
└─ README.md
