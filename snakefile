import glob
import os
import pandas as pd

configfile: "config.yaml"

IN_PATH = config["in_path"]
OUT_PATH = config["out_path"]
THREADS = config["threads"]
SAMPLES=config["samples"]
ITER = config["iterations"]

HAPS = [1, 2]

rule all:
    input:
        # expand(os.path.join(OUT_PATH,"02.Mapping/{sample}/{sample}.bam"), sample=SAMPLES),
        # expand(os.path.join(OUT_PATH,"03.Target/{sample}/{sample}_FMR1.bam"), sample=SAMPLES),
        # expand(os.path.join(OUT_PATH,"04.Assembly/{sample}/flye/assembly.fasta"), sample=SAMPLES),
        # expand(os.path.join(OUT_PATH,"05.Cor/{sample}/racon_round{round}.fasta"),round=[ITER], sample=SAMPLES),
        expand(os.path.join(OUT_PATH,"05.Cor/{sample}/quast_round{round}"),round=[ITER], sample=SAMPLES),
        # expand(os.path.join(OUT_PATH,"06.Alle/{sample}/{sample}.allele{hap}.bam"),sample=SAMPLES,hap=HAPS),
        os.path.join(OUT_PATH,"07.Summary/STR_summary.tsv"),


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

rule mapping:
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

rule extract_region:
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
        out_dir=directory(os.path.join(OUT_PATH,"04.Assembly/{sample}/flye")),
        out_contig=os.path.join(OUT_PATH,"04.Assembly/{sample}/flye/assembly.fasta"),
    threads: THREADS["high"]
    shell:
        """
        flye --nano-raw {input.r1} --genome-size 100k --threads {threads} --out-dir {output.out_dir}
        """

rule polish:
    input:
        r1=os.path.join(OUT_PATH,"01.QC/{sample}/{sample}.clean.fastq.gz"),
        ref=lambda w: os.path.join(OUT_PATH,f"05.Cor/{w.sample}/racon_round{int(w.round)-1}.fasta") if int(w.round) > 1 else os.path.join(OUT_PATH,f"04.Assembly/{w.sample}/flye/assembly.fasta"),
    output:
        paf=os.path.join(OUT_PATH,"05.Cor/{sample}/minimap2_round{round}.paf"),
        fasta=os.path.join(OUT_PATH,"05.Cor/{sample}/racon_round{round}.fasta"),
    threads: THREADS["high"]
    shell:
        """
        minimap2 -t {threads} -x map-ont {input.ref} {input.r1} > {output.paf}
        racon -t {threads} {input.r1} {output.paf} {input.ref} > {output.fasta}
        """

# ========== QUAST 评估 ==========
rule quast:
    input:
        fasta=lambda w: os.path.join(OUT_PATH,f"05.Cor/{w.sample}/racon_round{int(w.round)}.fasta")
    output:
        directory(os.path.join(OUT_PATH,"05.Cor/{sample}/quast_round{round}"))
    threads: THREADS["high"]
    shell:
        """
        quast {input.fasta} -o {output} -t {threads}
        """

rule re_mapping:
    input:
        # r1=os.path.join(OUT_PATH,"01.QC/{sample}/{sample}.clean.fastq.gz"),
        r1=os.path.join(OUT_PATH,"03.Target/{sample}/{sample}_FMR1.fq"),
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

rule haplotype:
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

rule re_polish:
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

rule str_loacte:
    input:
        allele=os.path.join(OUT_PATH,"06.Alle/{sample}/{sample}.allele{hap}.fa"),
        flank_file=config["flank_file"],
    output:
        bam=os.path.join(OUT_PATH,"06.Alle/{sample}/{sample}.allele{hap}.bam"),
    threads: THREADS["high"]
    params:
        chr=config["chr"]
    shell:
        r"""
        minimap2 -t {threads} -ax sr {input.allele} {input.flank_file} | samtools sort -o {output.bam}
        samtools index {output.bam}
        """

rule str_count:
    input:
        bam=os.path.join(OUT_PATH,"06.Alle/{sample}/{sample}.allele{hap}.bam"),
        allele=os.path.join(OUT_PATH,"06.Alle/{sample}/{sample}.allele{hap}.fa"),
    output:
        out_summary=os.path.join(OUT_PATH,"07.Summary/{sample}/{sample}_{hap}_STR_summary.tsv"),
    params:
        out_dir=directory(os.path.join(OUT_PATH,"07.Summary/{sample}")),
        no="{sample}_{hap}"
    shell:
        r"""
        python STR.py --bam {input.bam} --contig {input.allele} --sample {params.no} --out_dir {params.out_dir} 
        """

rule merge:
    input:
        expand(os.path.join(OUT_PATH,"07.Summary/{sample}/{sample}_{hap}_STR_summary.tsv"), sample=SAMPLES,hap=HAPS)
    output:
        os.path.join(OUT_PATH,"07.Summary/STR_summary.tsv")
    run:
        import pandas as pd
        dfs = [pd.read_csv(f, sep="\t") for f in input]
        df_all = pd.concat(dfs, ignore_index=True)
        df_all.to_csv(output[0], sep="\t", index=False)
