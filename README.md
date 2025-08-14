# Gene Annotator Pro - Offline High-Precision Annotation System

This tool is a Nextflow-based pipeline designed for comprehensive, offline functional annotation of transcriptome sequences. It integrates state-of-the-art bioinformatics tools and databases to produce systematic, scientifically rigorous, and informative annotation results.

## Key Features

-   **Fully Offline**: All analysis is performed without an internet connection once databases are prepared.
-   **Automated Preprocessing**: The pipeline automatically builds all necessary database indexes from user-provided raw files.
-   **Updatable Databases**: A simple command-line flag allows for the forced regeneration of all prepared databases from new raw files.
-   **Massively Parallel**: Designed to leverage multi-core servers and HPC clusters for maximum computational efficiency.
-   **Reproducible**: Utilizes Conda for strict software version management, ensuring scientific reproducibility.
-   **Intelligent Integration**: Employs a DuckDB backend for high-performance data integration, using a priority-based rule system to generate a final, human-readable report.

## Phase 1: Setup & Configuration

Before the first run, please complete the following setup.

### 1. Install Prerequisites

You must have `Conda` (or preferably `Mamba` for speed) and `Nextflow` installed.

```bash
# Install Mamba (recommended for faster environment creation)
conda install -n base -c conda-forge mamba

# Install Nextflow
mamba create -n nextflow -c bioconda nextflow

# Ensure you are in the 'nextflow' conda environment
conda activate nextflow

# Run the main pipeline
# The --input flag is required
# nextflow run main.nf -profile standard --input /path/to/your/transcripts.fasta --update_db
nextflow run main.nf -profile standard --input transcript.fasta 

# tree str
#.
#├── bin
#│   └── integrate_results.py
#├── database_prep.sh
#├── main.nf
#├── nextflow.config
#├── raw_database
#│   ├── eggnog_db
#│   ├── interproscan-5.75-106.0
#│   ├── interproscan-5.75-106.0-64-bit.tar.gz
#│   ├── interproscan-5.75-106.0-64-bit.tar.gz.md5
#│   ├── Pfam-A.clans.tsv.gz
#│   ├── Pfam-A.hmm.gz
#│   ├── Rfam.clanin
#│   ├── Rfam.cm
#│   ├── Rfam.cm.gz
#│   ├── Rfam.cm.i1f
#│   ├── Rfam.cm.i1i
#│   ├── Rfam.cm.i1m
#│   ├── Rfam.cm.i1p
#│   ├── uniprot_sprot.dmnd
#│   ├── uniprot_sprot.fasta.gz
#│   ├── uniprot_trembl.dmnd
#│   ├── uniprot_trembl.fasta.gz
#│   ├── uniref90.dmnd
#│   └── uniref90.fasta.gz
#├── README.md
#├── software_inst.sh
