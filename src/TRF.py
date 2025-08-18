import pandas as pd
import argparse

def trf_extract(input_file):
    # 筛选 period size = 3 的重复
    # 检查其 consensus motif 是否属于 CGG 的等价形式集合 {CGG, CCG, GGC, CGC, GCC, GCG}
    motifs = {"CGG", "CCG", "GGC", "CGC", "GCC", "GCG"}
    records = []
    with open(input_file, "r") as fin:
        for line in fin:
            if not line.strip():
                continue
            parts = line.strip().split()
            if not (parts[0].isdigit() and parts[1].isdigit()) or int(parts[2])%3 !=0 : # 判断前两列是否数字,第3列是否为3的倍数
                continue
            if parts[13] not in motifs: # # 判断第14列是否在motifs集合中
                continue
            records.append(parts)
    return records

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
    parser = argparse.ArgumentParser(description="CGG repeat counts of TRF software")
    parser.add_argument("-i", "--input", required=True, help="Input TRF bat file")
    parser.add_argument("-o","--output", required=True,help="out put file")
    parser.add_argument("--sample", required=True,help="sampleID")
    args = parser.parse_args() 
       
    records=trf_extract(args.input)    
    
    columns=[
        "start", "end", "period_size", "copy_number", "consensus_size", "percent_matches",
        "percent_indels", "alignment_score", "A", "C", "G", "T", "entropy", "motif", "consensus_sequence"
    ]
    if records:
        df = pd.DataFrame(records, columns=columns)
    else:
        df = pd.DataFrame(columns=columns)
    
    df["copy_number"] = df["copy_number"].astype(float)
    df["result"] = df["copy_number"].apply(classify)

    df.insert(0, "Sample", args.sample)
    col = df.pop('result')
    df.insert(1, 'result', col)

    df.to_csv(args.output, sep="\t", index=False)
    
if __name__ == "__main__":
    main()
