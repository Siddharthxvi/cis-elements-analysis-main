import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

def find_promoters(sequence):
    promoters = []
    search_start = 0
    seq_len = len(sequence)
    
    while True:
        n = sequence.find("TATA", search_start)
        if n == -1:
            break
        
        if n > 1000:
            found_atg = False
            for i in range(n + 28, n + 36):
                if i + 3 <= seq_len and sequence[i:i+3] == "ATG":
                    found_atg = True
                    break
            
            if found_atg:
                promoters.append(sequence[n - 1000: n])
        
        search_start = n + 1001
        if search_start >= seq_len:
            break
            
    return promoters

def calculate_gc(seq):
    seq_upper = seq.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    at_count = seq_upper.count('A') + seq_upper.count('T')
    total = gc_count + at_count
    if total > 0:
        return (gc_count / total) * 100.0, gc_count, at_count
    return 0.0, 0, 0

# Natural sorting for chromosomes (chr1, chr2, ... chr10)
def extract_numeric(filename):
    basename = os.path.basename(filename)
    num = ''.join(filter(str.isdigit, basename))
    return int(num) if num else float('inf')

def process_crop(crop_path):
    fasta_files = glob.glob(os.path.join(crop_path, "*.fasta"))
    if not fasta_files:
        return None
    
    fasta_files.sort(key=extract_numeric)
    
    genome_gc = 0
    genome_at = 0
    
    chrom_results = []
    
    for fasta_file in fasta_files:
        chrom_name = os.path.basename(fasta_file).replace('.fasta', '')
        
        with open(fasta_file, 'r') as f:
            sequence = []
            for line in f:
                line = line.strip()
                if not line.startswith('>'):
                    sequence.append(line.upper())
            full_seq = "".join(sequence)
            
            # Genome GC for this chr
            gc_pct, g_c, a_t = calculate_gc(full_seq)
            genome_gc += g_c
            genome_at += a_t
            
            # Promoters
            promoters = find_promoters(full_seq)
            if promoters:
                prom_seq = "".join(promoters)
                prom_gc_pct, _, _ = calculate_gc(prom_seq)
            else:
                prom_gc_pct = 0.0
                
            chrom_results.append({
                "Chromosome": chrom_name,
                "Promoter_GC_%": prom_gc_pct,
                "Num_Promoters": len(promoters)
            })
            
    total_genome_pct = 0.0
    if genome_gc + genome_at > 0:
        total_genome_pct = (genome_gc / (genome_gc + genome_at)) * 100.0
        
    return {
        "chrom_results": chrom_results,
        "genome_baseline": total_genome_pct
    }

def main():
    base_dir = "Data"
    if not os.path.exists(base_dir):
        print(f"Error: {base_dir} not found.")
        return
        
    crops = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    print(f"Found {len(crops)} crop directories in {base_dir}\n")
    
    all_tabular_data = []
    
    for crop in crops:
        crop_path = os.path.join(base_dir, crop)
        print(f"Processing {crop}...")
        
        crop_data = process_crop(crop_path)
        if crop_data is None:
            print(f"  -> Skipping {crop}, no fasta files found.")
            continue
            
        genome_baseline = crop_data["genome_baseline"]
        chrom_results = crop_data["chrom_results"]
        
        if not chrom_results:
            continue
            
        # Add to tabular data
        for cr in chrom_results:
            all_tabular_data.append({
                "Crop": crop,
                "Chromosome": cr["Chromosome"],
                "Promoter_GC_%": round(cr["Promoter_GC_%"], 2),
                "Num_Promoters": cr["Num_Promoters"],
                "Genome_Baseline_GC_%": round(genome_baseline, 2)
            })
            
        # Plotting
        df_chroms = pd.DataFrame(chrom_results)
        
        plt.figure(figsize=(10, 6))
        plt.bar(df_chroms["Chromosome"], df_chroms["Promoter_GC_%"], color='skyblue', label='Promoter GC %')
        plt.axhline(y=genome_baseline, color='red', linestyle='--', label=f'Genome Baseline ({genome_baseline:.2f}%)')
        
        plt.xlabel('Chromosome')
        plt.ylabel('GC Content (%)')
        plt.title(f'Promoter GC Content per Chromosome - {crop}')
        plt.xticks(rotation=45, ha='right')
        
        # fixed y limit for consistent comparison across crops
        plt.ylim(0, 50)
        plt.legend()
        plt.tight_layout()
        
        plot_filename = f"GC_Content_{crop}.png"
        plt.savefig(plot_filename)
        plt.close()
        print(f"  -> Saved plot to {plot_filename}")
        
    # Tabular output
    df_all = pd.DataFrame(all_tabular_data)
    df_all.to_csv("GC_Content_Summary.csv", index=False)
    print("\nSaved all tabular data to GC_Content_Summary.csv")

if __name__ == "__main__":
    main()
