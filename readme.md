# FMR1_STR pipeline

- update at     2025-08-18 
- author        Hu Qiaoli  

## 1. Introduction
The **FMR1_STR pipeline** is designed to analyze **FMR1 gene 5’UTR CGG repeat expansions** associated with **Fragile X Syndrome** using **Nanopore Technology long reads**.  

The pipeline integrates **quality control, mapping, target region assembly, haplotype-specific phasing, and STR repeat counting**.  This allows for accurate detection of repeat length, motif interruptions, and allele-specific haplotypes.  

**Pipeline repository:**
```bash
git clone https://github.com/<YourRepo>/FMR1_STR.git
```

## 2. Requirements
**Docker**

We recommend running this pipeline inside a docker or Apptainer container.
Our docker image contains the following bioinformatics tools :

1. NanoPlot (QC visualization)
2. bedtools (region extraction)
3. filtlong (read filtering )
4. flye (de novo assembly of STR region)
5. longshot (SNV-based haplotype phasing)
6. minimap2 (long-read mapping)
7. quast (assembly QC)
8. racon (contig polishing)
9. samtools (alignment handling)
10. snakemake (pipeline execution engine)
11. Python3 (v3.10.16) is bundled with the following packages:
    1. numpy
    2. pandas
    3. argparse

## 3. Script Structure
```
./FMR1_STR/
├── config.yaml                 ## snakemake config file
├── Dockerfile                  ## Dockerfile for docker image
├── lib/                        ## 
│   ├── STR.py                  ## STR repeat annotation script
│   ├── hg38_chrXY.*            ## reference  (hg38 chrX subset)
│   ├── fmr1_flanks.fa          ## FMR1 STR flanking sequences
├── readme.txt                  ## documentation
└── run.py                      ## main entry of the 
├── readme.md
├── snakefile                   ## snakemake file
└── workflow.pngs
```

## 4. Usage
Run the pipeline via the provided shell script (example):
```
snakemake -np                   ## dry-run 
snakemake -j ${core_num} -p 
```
docker run -idtv /gpfsvol1/P_BIOINFO/USER/huql/:/gpfsvol1/P_BIOINFO/USER/huql/ --name=${name} registry.servicemgr.gendow:5000/huql/meta:v2 /bin/bash


Input Files
-in : FASTQ directory containing raw nanopore reads


## 5. Output Directory
```
./output/
├── 01.QC/             ### Raw and filtered reads QC
├── 02.Mapping/        ### Reads mapped to reference
├── 03.Target/         ### Reads aligned to FMR1 region
├── 04.Assembly/       ### De novo assembly of STR region
├── 05.Cor/            ### Polished contig
├── 06.Alle/           ### Haplotype contigs after phasing
├── FMR1.bed           ### STR locus BED file
└── report.txt         ### STR size, motif, and summary statistics
```

## 6. Notes
- For reproducibility, always run inside Docker/Apptainer.

- The pipeline is modular: each step can be re-run independently.

- Downstream manual inspection (IGV + TRF reports) is recommended for final confirmation.

- Multi-sample runs can be managed by Snakemake batch mode.

## 7. Workflow
![alt text](workflow.png)