import pysam
from Bio import SeqIO
import sys
import os
import argparse
import pandas as pd

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def get_contig_seq(fasta_file, seq_name):
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id == seq_name:
            return str(record.seq).upper()
    return None
     
def get_str_region(bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    coords = {"left_flank": None, "right_flank": None}

    for aln in bam.fetch():
        qname = aln.query_name
        if qname not in coords:
            continue
        if not aln.is_unmapped:
            coords[qname] = (aln.reference_name, 
                             aln.reference_start, 
                             aln.reference_end, 
                             "-" if aln.is_reverse else "+")

    bam.close()
    # 如果都没比对上
    if coords["left_flank"] is None and coords["right_flank"] is None:
        return None
    return coords

def get_str_seq(coords,contig_fasta):
    if coords is None:
        # print("No flank sequences found in the BAM file.")
        return None, None, None, None
        # sys.exit(0)
    left_flank = coords["left_flank"]
    right_flank = coords["right_flank"]
    
    if left_flank and right_flank and left_flank[0]==right_flank[0]:
        str_start = left_flank[2]  # 左侧翼末端
        str_end = right_flank[1]   # 右侧翼起点
        contig_name=left_flank[0]
        str_seq=get_contig_seq(contig_fasta, contig_name)[str_start:str_end]
    elif left_flank and right_flank and (left_flank[0]!=right_flank[0]):
        print("Warning: Left and right flanks map to different contigs.")
        contig_name=left_flank[0]
        contig_seq=get_contig_seq(contig_fasta, contig_name)
        str_start = left_flank[2] # 取左侧翼所在的contig，默认1000bp窗口
        str_end = min(len(contig_seq), str_start+1000)
        str_seq=contig_seq[str_start:str_end]
    elif left_flank and not right_flank:
        contig_name=left_flank[0]
        contig_seq=get_contig_seq(contig_fasta, contig_name)
        str_start = left_flank[2]
        str_end = min(len(contig_seq), str_start+1000)  # 默认1000bp窗口
        str_seq=contig_seq[str_start:str_end]
    elif not left_flank and right_flank:
        contig_name=right_flank[0]
        contig_seq=get_contig_seq(contig_fasta, contig_name)
        str_end = right_flank[1]
        str_start = max(0, str_end-1000)
        str_seq=contig_seq[str_start:str_end]
    else:
        return None, None, None, None
    
    return contig_name,str_start,str_end,str_seq

# -----------------------
# CGG repeat检测函数（允许AGG/CAG中断、1-2 bp mismatch）
# -----------------------
def count_cgg_repeat(seq, max_mismatch=1, allow_agg=True):
    if seq is None or len(seq) < 3:
        return 0, ""
    i = 0
    repeat_count = 0
    motif_seq = []
    while i <= len(seq)-3:
        tri = seq[i:i+3]
        # 判断是否CGG或可接受的中断
        mismatch = sum(1 for a,b in zip(tri,"CGG") if a!=b)
        if mismatch <= max_mismatch or (allow_agg and tri in ["AGG","CAG"]):
            repeat_count +=1
            motif_seq.append(tri)
            i += 3
        else:
            i +=1             # 连续CGG断裂就跳到下一个碱基
    return repeat_count, "".join(motif_seq)

def classify(count):
    if count < 45:
        return "Normal"
    elif 45 <= count <= 54:
        return "Intermediate"
    elif 55 <= count <= 200:
        return "Premutation"
    else:
        return "FullMutation"

def main():
    parser = argparse.ArgumentParser(description="Scan BAM for CGG repeat counts")
    parser.add_argument("-b", "--bam", required=True, help="Input BAM file")
    parser.add_argument("--contig", required=True, help="contig fasta file")
    parser.add_argument("--max_mm", type=int, default=1)
    parser.add_argument("--out_dir", required=True,help="out put directory")
    parser.add_argument("--sample", required=True,help="sampleID")
    args = parser.parse_args()
    
    os.makedirs(args.out_dir, exist_ok=True)
    
    coords=get_str_region(args.bam)
    contig_name,str_start,str_end,str_seq=get_str_seq(coords, args.contig)
    
    if str_seq is not None:
        repeat_count, str_motifs = count_cgg_repeat(str_seq, max_mismatch=args.max_mm, allow_agg=True) # True 表示AGG/CAG允许中断
        rev_repeat_count, rev_str_motifs = count_cgg_repeat(revcomp(str_seq), max_mismatch=args.max_mm, allow_agg=True) # True 表示AGG/CAG允许中断
        result=[classify(repeat_count),classify(rev_repeat_count)]
    else:
        repeat_count = rev_repeat_count = "NA"
        str_motifs = rev_str_motifs = ""
        result = ["NA", "NA"]

    df = pd.DataFrame([[
        args.sample,
        contig_name, 
        str_start, str_end, 
        len(str_seq) if str_seq else 0 , 
        result,
        repeat_count, 
        rev_repeat_count, 
        str_motifs,
        rev_str_motifs, 
        str_seq]],
        columns=[
            "sampleID",
            "contig",
            "STR_start",
            "STR_end",
            "str_len",
            "result",
            "repeat_count",
            "rev_repeat_count", 
            "motifs",
            "rev_motifs",
            "STR_sequence"
            ])
    df.to_csv(f"{args.out_dir}/{args.sample}_STR_summary.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
