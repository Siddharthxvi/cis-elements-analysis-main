import os
import glob
import pandas as pd
import numpy as np

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
        
        search_b_from = l + 1 if a == b else l
        h = promoter_sequence.find(b, search_b_from)
        
        if h == -1: break
        
        a_end = l + len(a)
        raw_spacer_length = h - a_end
        
        # Use spacer 0 to signify adjacent or overlapping motif pairs
        if raw_spacer_length < 0:
            reported_len = 0
        else:
            reported_len = raw_spacer_length
            
        spacers.append({
            'len': reported_len,
            'is_overlapping': raw_spacer_length < 0,
            'is_adjacent': raw_spacer_length == 0
        })
        
        # move pos to next potential valid site
        pos = h + len(b)
        if pos <= l: pos = l + 1
        
    return spacers

def calculate_stats(lengths):
    if not lengths:
        return np.nan, np.nan, np.nan, np.nan, np.nan
        
    s_mean = np.mean(lengths)
    s_median = np.median(lengths)
    s_sd = np.std(lengths, ddof=1) if len(lengths) > 1 else 0.0
    
    q1 = np.percentile(lengths, 25)
    q3 = np.percentile(lengths, 75)
    s_iqr = q3 - q1
    
    # Mode
    vals, counts = np.unique(lengths, return_counts=True)
    max_count = np.max(counts)
    modes = vals[counts == max_count]
    s_mode = np.min(modes) # select the smallest mode dynamically
    
    return s_mean, s_median, s_mode, s_sd, s_iqr

def main():
    print("=== DistroStats Computations ===")
    base_dir = "Data"
    if not os.path.exists(base_dir): 
        print("Data directory not found.")
        return
        
    crops = sorted([d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))])
    motif_pairs = [("AAAG", "AAAG"), ("AAAG", "ACGT"), ("ACGT", "AAAG"), ("ACGT", "ACGT")]
    
    results = []
    
    for crop in crops:
        crop_path = os.path.join(base_dir, crop)
        fasta_files = sorted(glob.glob(os.path.join(crop_path, "*.fasta")))
        if not fasta_files: continue
        
        # Extract promoters aggressively
        sequences = []
        for f in fasta_files:
            with open(f, 'r') as file:
                seq_lines = []
                for line in file:
                    if not line.startswith('>'):
                        seq_lines.append(line.strip().upper())
                sequences.append("".join(seq_lines))
                
        promoters = []
        for seq in sequences:
            promoters.extend(find_promoters(seq))

        for element_a, element_b in motif_pairs:
            print(f"Processing -> Crop: {crop} | Motifs: {element_a}-{element_b}")
            
            all_spacers = []
            for prom in promoters:
                all_spacers.extend(find_spacers(prom, element_a, element_b))
                
            n_instances = len(all_spacers)
            if n_instances == 0:
                continue
                
            lengths = [s['len'] for s in all_spacers]
            mean_s, med_s, mode_s, sd_s, iqr_s = calculate_stats(lengths)
            
            zero_count = sum(1 for s in all_spacers if s['len'] == 0)
            zero_frac = zero_count / n_instances
            has_composite = zero_count > 0
            
            overlapping_count = sum(1 for s in all_spacers if s['is_overlapping'])
            adjacent_count = sum(1 for s in all_spacers if s['is_adjacent'])
            
            filtered_lengths = [l for l in lengths if l > 0]
            mean_nz, med_nz, _, sd_nz, iqr_nz = calculate_stats(filtered_lengths)
            
            results.append({
                "Crop": crop,
                "CRE_type": f"{element_a}-{element_b}",
                "N_instances": n_instances,
                "Mean_spacer": round(mean_s, 2),
                "Median_spacer": med_s,
                "Mode_spacer": f"{mode_s} (composite driven mode!)" if mode_s == 0 else mode_s,
                "SD_spacer": round(sd_s, 2),
                "IQR_spacer": round(iqr_s, 2),
                "Zero_spacer_count": zero_count,
                "Zero_spacer_fraction": round(zero_frac, 4),
                "Has_composite_elements": has_composite,
                "Overlapping_count": overlapping_count,
                "Adjacent_count": adjacent_count,
                "Mean_spacer_no_zero": round(mean_nz, 2) if not np.isnan(mean_nz) else "N/A",
                "Median_spacer_no_zero": med_nz if not np.isnan(med_nz) else "N/A",
                "SD_spacer_no_zero": round(sd_nz, 2) if not np.isnan(sd_nz) else "N/A",
                "IQR_spacer_no_zero": round(iqr_nz, 2) if not np.isnan(iqr_nz) else "N/A"
            })
            
    df = pd.DataFrame(results)
    output_filename = "Spacer_Distribution_Statistics.csv"
    df.to_csv(output_filename, index=False)
    
    print(f"\n\nAnalysis Complete -> Saved robust table to {output_filename}")

if __name__ == "__main__":
    main()
