# FMR1_STR pipeline

- update at     2025-08-18 
- author        Hu Qiaoli 
- version       0.1

## 1. Introduction
The **FMR1_STR pipeline** is designed to analyze **FMR1 gene 5’UTR CGG repeat expansions** which are associated with **Fragile X Syndrome(FXS)** using **Nanopore Technology long reads**.  

The pipeline integrates **quality control, mapping, target region extraction and de novo assembly, phasing of haplotype, STR repeat counting and Summary report generation**.  This allows for accurate detection of repeat length, motif interruptions, and allele-specific haplotypes.  

**Pipeline repository:**
```bash
git clone https://github.com/ivyhulife/FMR1_STR.git
```

## 2. Requirements

We recommend running this pipeline inside a **docker** or **Apptainer** container for full reproducibility. The pre-built Docker image includes the following tools:
| Tool          | Version      | Purpose                                    |
| ------------- | ------------ | ------------------------------------------ |
| **NanoPlot**  | 1.32.0       | QC visualization of Nanopore reads         |
| **bedtools**  | 2.31.1      | Read extraction from genomic regions       |
| **filtlong**  | 0.2.1       | Filtering of Nanopore reads                |
| **flye**      | 2.9.3-b1797 | De novo assembly of STR region             |
| **longshot**  | 1.0.0       | SNV-based haplotype phasing                |
| **minimap2**  | 2.30-r1287    | Long-read mapping                          |
| **quast**     | 5.3.0         | Assembly quality assessment                |
| **racon**     | 1.5.0       | Contig polishing                           |
| **samtools**  | 1.21        | Alignment manipulation                     |
| **snakemake** | 7.32.4        | Workflow management                        |
| **Python3**   | 3.10.16      | Workflow scripts (numpy, pandas, argparse,pysam,biopython) |

## 3. Script Directory Structure
```
FMR1_STR/
├── config.yaml                 # Global configuration for Snakemake
├── demo_results.png            # The results of demo reads
├── Dockerfile                  # Docker image definition
├── lib/                        
│   ├── fmr1_flanks.fa          # Flanking sequences for STR localization
│   ├── hg38_chrXY.*            # Reference (hg38 chromosome X subset)
│   ├── STR.py                  # STR repeat detection and annotation
├── readme.md                   # Project documentation
├── snakefile                   # Snakemake pipeline definition
└── workflow.png                # Workflow schematic
```

## 4. Usage
### 4.1. Edit **config.yaml** before running the pipeline:
```yaml
samples:  ["SRR17138637","SRR17138639"]   # List of sample IDs
in_path: "./raw_data/"           # Directory containing input FASTQ files
out_path: "./results/"             # Output directory     

threads:
  normal: 16                      #  recommended for small or lightweight steps
  high: 60                        # recommended for computationally intensive step

```
### 4.2. Running the Pipeline via snakemeke :
```
snakemake -np                   ## Dry-run (no execution, check DAG)
snakemake -j ${core_num} -p     ## Full execution

# Run inside Docker
docker build -t fmr1_str:latest .
docker run -itv $(pwd):/data fmr1_str:latest  snakemake -np 
docker run -itv $(pwd):/data fmr1_str:latest  snakemake -j ${core_num} -p 
```
### 4.3. Workflow
![alt text](workflow.png)

### 4.4. STR_summary report in ./output/07.Summary
![alt text](demo_results.png)

## 5. Output Directory Structure
```
./output/
├── 01.QC/             # QC reports for raw and filtered reads
├── 02.Mapping/        # Reads mapped to reference (hg38 subset)
├── 03.Target/         # Extracted reads covering the FMR1 locus
├── 04.Assembly/       # De novo assembly of STR region
├── 05.Cor/            # Polished contig
├── 06.Alle/           # Haplotype-specific contigs after phasing
├── 07.Summary/        # STR size, motif structure, summary statistics
└── FMR1.bed           # STR locus coordinates
```

## 6. Notes
- Always run the pipeline inside Docker/Apptainer for reproducibility.

- The workflow is modular: each step can be executed independently by specifying the rule name in Snakemake.

- Multi-sample runs can be managed by Snakemake batch mode.