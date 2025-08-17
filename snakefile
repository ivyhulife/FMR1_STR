import glob
import os
import pandas as pd

configfile: "config.yaml"
IN_PATH = config["in_path"]
OUT_PATH = config["out_path"]
THREADS = config["threads"]
SAMPLES=config["samples"]
ITER = config["iterations"]
# SAMPLES= [os.path.basename(f).replace(".fastq.gz","") for f in glob.glob(os.path.join(IN_PATH,"*.fastq.gz"))]
HAPS = [1, 2]

rule all:
    input:
        # expand(os.path.join(OUT_PATH,"02.Mapping/{sample}/{sample}.bam"), sample=SAMPLES),
        # expand(os.path.join(OUT_PATH,"03.Target/{sample}/{sample}_FMR1.bam"), sample=SAMPLES),
        expand(os.path.join(OUT_PATH,"04.assembly/{sample}/flye/assembly.fasta"), sample=SAMPLES),
        expand(os.path.join(OUT_PATH,"05.Cor/{sample}/racon_round{round}.fasta"),round=[ITER], sample=SAMPLES),
        expand(os.path.join(OUT_PATH,"06.Alle/{sample}/{sample}.allele{hap}.fa"),sample=SAMPLES,hap=HAPS),

rule create_bed:
    output:
        bed_file = os.path.join(OUT_PATH,"FMR1.bed")
    run:
        chr=config["chr"]
        start = config["fmr1_start"] - config.get("pad_bp",50000)
        end = config["fmr1_start"] + config.get("pad_bp",50000)
        shell(f"echo -e '{chr}\t{start}\t{end}' > {output.bed_file}")

rule qc:
    input:
        r1=lambda wildcards: os.path.join(IN_PATH, f"{wildcards.sample}.fastq.gz"),
    output:
        raw_qc=directory(os.path.join(OUT_PATH,"01.QC/{sample}/raw_nanoplot")),
        r1 = os.path.join(OUT_PATH,"01.QC/{sample}/{sample}.clean.fastq.gz"),
        clean_qc=directory(os.path.join(OUT_PATH,"01.QC/{sample}/after_nanoplot")),
    threads: THREADS["high"]
    shell:
        """
        /ifs1/Software/miniconda3/envs/nanoplot/bin/NanoPlot -t {threads} --fastq {input.r1} -o {output.raw_qc}
        filtlong --min_length 100 --min_mean_q 80 {input.r1} |  gzip > {output.r1} 
        /ifs1/Software/miniconda3/envs/nanoplot/bin/NanoPlot -t {threads} --fastq {output.r1} -o {output.clean_qc}
        """

rule minimap2:
    input:
        r1=os.path.join(OUT_PATH,"01.QC/{sample}/{sample}.clean.fastq.gz"),
        ref_file=config["ref_path"]
    output:
        bam=os.path.join(OUT_PATH,"02.Mapping/{sample}/{sample}.bam"),
    threads: THREADS["high"]
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.ref_file} {input.r1} | samtools sort -o {output.bam}
        samtools index -@ {threads} {output.bam}
        """

rule target:
    input:
        bam=os.path.join(OUT_PATH,"02.Mapping/{sample}/{sample}.bam"),
        bed_file = os.path.join(OUT_PATH,"FMR1.bed")
    output:
        bam=os.path.join(OUT_PATH,"03.Target/{sample}/{sample}_FMR1.bam"),
        depth=os.path.join(OUT_PATH,"03.Target/{sample}/{sample}_FMR1.depth.txt"),
        cov=os.path.join(OUT_PATH,"03.Target/{sample}/{sample}_FMR1.cov.txt"),
        fq=os.path.join(OUT_PATH,"03.Target/{sample}/{sample}_FMR1.fq"),
    threads: THREADS["high"]
    shell:
        """
        samtools view -@ {threads} -b -L {input.bed_file} {input.bam} > {output.bam}
        samtools index -@ {threads} {output.bam}
        samtools depth -@ {threads} -b {input.bed_file} {output.bam} > {output.depth}
        bedtools coverage -a {input.bed_file} -b {output.bam} > {output.cov}
        samtools fastq -@ {threads} {output.bam} > {output.fq}
        """

rule assembly:
    input:
        r1=os.path.join(OUT_PATH,"03.Target/{sample}/{sample}_FMR1.fq"),
    output:
        out_dir=directory(os.path.join(OUT_PATH,"04.assembly/{sample}/flye")),
        out_contig=os.path.join(OUT_PATH,"04.assembly/{sample}/flye/assembly.fasta"),
    threads: THREADS["high"]
    shell:
        """
        flye --nano-raw {input.r1} --genome-size 100k --threads {threads} --out-dir {output.out_dir}
        """

rule polish1:
    input:
        r1=os.path.join(OUT_PATH,"01.QC/{sample}/{sample}.clean.fastq.gz"),
        ref=lambda w: os.path.join(OUT_PATH,f"05.Cor/{w.sample}/racon_round{int(w.round)-1}.fasta") if int(w.round) > 1 else os.path.join(OUT_PATH,f"04.assembly/{w.sample}/flye/assembly.fasta"),
    output:
        paf=os.path.join(OUT_PATH,"05.Cor/{sample}/minimap2_round{round}.paf"),
        fasta=os.path.join(OUT_PATH,"05.Cor/{sample}/racon_round{round}.fasta"),
    threads: THREADS["high"]
    shell:
        """
        minimap2 -t {threads} -x map-ont {input.ref} {input.r1} > {output.paf}
        racon -t {threads} {input.r1} {output.paf} {input.ref} > {output.fasta}
        """

# # ========== MEDAKA ==========
# rule medaka:
#     input:
#         r1=os.path.join(OUT_PATH,"01.QC/{sample}/{sample}.clean.fastq.gz"),
#         ref=os.path.join(OUT_PATH,"05.Cor/{sample}/racon_round{ITER}.fasta")
#     output:
#         out_dir=directory(os.path.join(OUT_PATH,"05.Cor/{sample}"))
#         consensus=os.path.join(OUT_PATH,"05.Cor/{sample}/{sample}_consensus.fasta")
#     threads: THREADS["high"]
#     shell:
#         """
#         medaka_consensus -i {input.r1} -d {input.ref} \
#             -o {output.out_dir} -t {threads} -m {MEDAKA_MODEL}
#         """

# # ========== QUAST 评估 ==========
# rule quast:
#     input:
#         fasta=lambda w: os.path.join(OUT_PATH,f"05.Cor/{w.sample}/racon_round{int(w.round)}.fasta")
#     output:
#         directory(os.path.join(OUT_PATH,"05.Cor/{sample}/quast_round{round}"))
#     threads: THREADS["high"]
#     shell:
#         """
#         quast {input.fasta} -o {output} -t {threads}
#         """

rule minimap2_target:
    input:
        r1=os.path.join(OUT_PATH,"01.QC/{sample}/{sample}.clean.fastq.gz"),
        contig=lambda w:os.path.join(OUT_PATH,f"05.Cor/{w.sample}/racon_round{ITER}.fasta")
    output:
        bam=os.path.join(OUT_PATH,"06.Alle/{sample}/{sample}.bam"),
    threads: THREADS["high"]
    shell:
        """
        minimap2 -t {threads} -ax map-ont {input.contig} {input.r1} | samtools sort -o {output.bam}
        samtools index -@ {threads} {output.bam}
        samtools faidx {input.contig}
        """

rule haplotype_assembly:
    input:
        contig=lambda w:os.path.join(OUT_PATH,f"05.Cor/{w.sample}/racon_round{ITER}.fasta"),
        bam=os.path.join(OUT_PATH,"06.Alle/{sample}/{sample}.bam"),
    output:
        vcf=os.path.join(OUT_PATH,"06.Alle/{sample}/{sample}.vcf"),
        hp_bam=os.path.join(OUT_PATH,"06.Alle/{sample}/{sample}.hp.bam"),
    threads: THREADS["high"]
    shell:
        r"""
        longshot --bam {input.bam} --ref {input.contig} --out {output.vcf} --out_bam {output.hp_bam}
        """

rule polish2:
    input:
        contig=lambda w:os.path.join(OUT_PATH,f"05.Cor/{w.sample}/racon_round{ITER}.fasta"),
        hp_bam=os.path.join(OUT_PATH,"06.Alle/{sample}/{sample}.hp.bam"),
    output:
        fq=os.path.join(OUT_PATH,"06.Alle/{sample}/{sample}.hp{hap}.fq"),
        paf=os.path.join(OUT_PATH,"06.Alle/{sample}/{sample}.hp{hap}.paf"),
        allele=os.path.join(OUT_PATH,"06.Alle/{sample}/{sample}.allele{hap}.fa"),
    threads: THREADS["high"]
    shell:
        r"""
        samtools view -bh -d HP:{wildcards.hap} {input.hp_bam} | samtools fastq - > {output.fq}
        minimap2 -x map-ont -t {threads} {input.contig} {output.fq} > {output.paf}
        racon -t {threads} {output.fq} {output.paf} {input.contig} > {output.allele}
        """

# ========== 修改中 ==========
rule extract:
    input:
        r1=lambda w: os.path.join(OUT_PATH,f"05.Cor/{w.sample}/racon_round{int(w.round)}.fasta"),
        ref_file=config["ref_path"],
        bed_file = os.path.join(OUT_PATH,"FMR1.bed"),
    output:
        bam=os.path.join(OUT_PATH,"05.Cor/{sample}/map_round{round}.bam"),
        tar_bam=os.path.join(OUT_PATH,"05.Cor/{sample}/map_round{round}_FMR1.bam"),
    threads: THREADS["high"]
    params:
        chr=config["chr"]
    shell:
        r"""
        set -euo pipefail
        minimap2 -t {threads} -ax asm20 {input.ref_file} {input.r1} | samtools sort -o {output.bam}
        samtools index {output.bam}
        samtools view -@ {threads} -b -L {input.bed_file} {output.bam} > {output.tar_bam}
        samtools index {output.tar_bam}
        """

rule target_fq:
    input:
        r1=lambda w: os.path.join(OUT_PATH,f"05.Cor/{w.sample}/racon_round{int(w.round)}.fasta"),
        tar_bam=lambda w:os.path.join(OUT_PATH,f"05.Cor/{w.sample}/map_round{int(w.round)}_FMR1.bam"),
    output:
        fasta=os.path.join(OUT_PATH,"05.Cor/{sample}/round{round}_FMR1.fasta"),
    threads: THREADS["high"]
    params:
        chr=config["chr"]
    shell:
        r"""
        min_start=$(samtools view -@ {threads} {input.tar_bam} | cut -f 6 | cut -d 'H' -f1 | sort -n | head -n1)
        start=$((min_start - 50000))
        end=$((min_start + 90000))
        if [ $start -lt 1 ]; then start=1; fi
        seqkit faidx {input.r1} {params.chr}:${{start}}-${{end}} | seqtk seq -A > {output.fasta}
        """

rule STR_count:
    input:
        bam=os.path.join(OUT_PATH,"03.Target/{sample}/{sample}_FMR1.bam"),
    output:
        out_dir=directory(os.path.join(OUT_PATH,"04.STR/{sample}")),
        # hist=f"{outdir}/{sample}.FMR1_CGG_hist.png",
        # report=f"{outdir}/{sample}_FMR1_CGG_report.txt"
    params:
        chr=config["chr"],
        start = config["fmr1_start"] - config.get("pad_bp",50000),
        end = config["fmr1_start"] + config.get("pad_bp",50000),
    shell:
        """
        python STR.py --bam {input.bam} \
            --chr {params.chr} --start {params.start} --end {params.end} \
            --out_dir {output.out_dir}
        """
