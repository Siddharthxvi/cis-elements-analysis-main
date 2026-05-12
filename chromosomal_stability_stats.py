import os
import glob
import numpy as np
import pandas as pd
from scipy import stats

def find_promoters(sequence):
    promoters = []
    search_start = 0
    seq_len = len(sequence)
    while True:
        n = sequence.find("TATA", search_start)
        if n == -1: break
        if n > 1000:
            found_atg = False
            for i in range(n + 28, n + 36):
                if i + 3 <= seq_len and sequence[i:i+3] == "ATG":
                    found_atg = True
                    break
            if found_atg:
                promoters.append(sequence[n - 1000: n])
        search_start = n + 1001
        if search_start >= seq_len: break
    return promoters

def find_spacers(promoter_sequence, a, b):
    spacers = []
    n = len(promoter_sequence)
    pos = 0
    while pos < n:
        l = promoter_sequence.find(a, pos)
        if l == -1: break
        search_start = l + len(a)
        if search_start >= n: break
        h = promoter_sequence.find(b, search_start)
        if h == -1: break
        space_len = h - search_start
        spacers.append(promoter_sequence[search_start:search_start+space_len])
        pos = h + len(b)
        if pos <= l: pos = l + 1
    return spacers

def calculate_gc(seq):
    seq_upper = seq.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    at_count = seq_upper.count('A') + seq_upper.count('T')
    total = gc_count + at_count
    if total > 0: return (gc_count / total) * 100.0
    return float('nan')

def main():
    print("=== Chromosomal Stability Statistics ===")
    
    # We will test the specific motif pairs commonly requested
    motif_pairs = [("AAAG", "AAAG"), ("AAAG", "ACGT"), ("ACGT", "AAAG"), ("ACGT", "ACGT")]
    base_dir = "Data"
    
    if not os.path.exists(base_dir):
        print(f"Error: {base_dir} not found.")
        return
        
    crops = sorted([d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))])
    
    results = []

    for element_a, element_b in motif_pairs:
        print(f"\\nAnalyzing {element_a}-{element_b} ...")
        
        for crop in crops:
            crop_path = os.path.join(base_dir, crop)
            fasta_files = sorted(glob.glob(os.path.join(crop_path, "*.fasta")))
            if not fasta_files: continue
            
            # To store individual spacer GC percentages keyed by chromosome
            chrom_spacer_gcs = {}
            # To store just the single mean GC% representing the whole chromosome
            chrom_means = []
            
            for fasta_file in fasta_files:
                chrom_name = os.path.basename(fasta_file).replace('.fasta', '')
                
                with open(fasta_file, 'r') as f:
                    sequence = []
                    for line in f:
                        line = line.strip()
                        if not line.startswith('>'):
                            sequence.append(line.upper())
                    full_seq = "".join(sequence)
                    
                    promoters = find_promoters(full_seq)
                    
                    # Store all individual spacers for this chromosome
                    current_chrom_gc_list = []
                    
                    total_chrom_g_c = 0
                    total_chrom_a_t = 0
                    
                    for prom in promoters:
                        spacers = find_spacers(prom, element_a, element_b)
                        for spacer in spacers:
                            seq_upper = spacer.upper()
                            g_c = seq_upper.count('G') + seq_upper.count('C')
                            a_t = seq_upper.count('A') + seq_upper.count('T')
                            if (g_c + a_t) > 0:
                                gc_pct = (g_c / (g_c + a_t)) * 100.0
                                current_chrom_gc_list.append(gc_pct)
                                total_chrom_g_c += g_c
                                total_chrom_a_t += a_t
                    
                    if current_chrom_gc_list:
                        chrom_spacer_gcs[chrom_name] = current_chrom_gc_list
                        chrom_mean_gc = (total_chrom_g_c / (total_chrom_g_c + total_chrom_a_t)) * 100.0
                        chrom_means.append(chrom_mean_gc)

            if not chrom_means:
                continue
                
            # 1. Coefficient of Variation across chromosomes
            # Formula: (Standard Deviation of Chromosome Means / Mean of Chromosome Means) * 100
            if len(chrom_means) > 1:
                cv = (np.std(chrom_means, ddof=1) / np.mean(chrom_means)) * 100
            else:
                cv = 0.0
                
            # 2. One-way ANOVA 
            # Null hypothesis: The means of individual spacer GC contents are equal across all chromosomes
            anova_f = float('nan')
            anova_p = float('nan')
            
            # Prepare groups (lists of GC% values per chromosome)
            groups = list(chrom_spacer_gcs.values())
            # Ensure we have at least 2 chromosomes and they have variability to perform ANOVA
            if len(groups) > 1 and sum(len(g) for g in groups) > len(groups):
                try:
                    f_stat, p_val = stats.f_oneway(*groups)
                    anova_f = f_stat
                    anova_p = p_val
                except Exception as e:
                    pass

            results.append({
                "Motif_Pair": f"{element_a}-{element_b}",
                "Crop": crop,
                "Num_Chromosomes": len(chrom_means),
                "Overall_Mean_GC_%": round(np.mean(chrom_means), 2),
                "CV_Across_Chromosomes_%": round(cv, 4),
                "ANOVA_F_Statistic": round(anova_f, 4) if not np.isnan(anova_f) else "N/A",
                "ANOVA_p_value": round(anova_p, 6) if not np.isnan(anova_p) else "N/A"
            })
            
    df = pd.DataFrame(results)
    csv_file = "Chromosomal_Stability_Analysis.csv"
    df.to_csv(csv_file, index=False)
    
    print("\n\n=== Final Results Summary ===")
    print(df.to_string())
    print(f"\nSaved statistical summary to {csv_file}")

if __name__ == "__main__":
    main()
