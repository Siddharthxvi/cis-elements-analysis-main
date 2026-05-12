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

def find_spacers(promoter_sequence, a, b):
    spacers = []
    n = len(promoter_sequence)
    pos = 0
    
    while pos < n:
        l = promoter_sequence.find(a, pos)
        if l == -1:
            break
            
        search_start = l + len(a)
        if search_start >= n:
            break
            
        h = promoter_sequence.find(b, search_start)
        if h == -1:
            break
            
        space_len = h - search_start
        spacers.append(promoter_sequence[search_start:search_start+space_len])
        
        pos = h + len(b)
        if pos <= l: 
            pos = l + 1
            
    return spacers

def calculate_gc(seq):
    seq_upper = seq.upper()
    gc_count = seq_upper.count('G') + seq_upper.count('C')
    at_count = seq_upper.count('A') + seq_upper.count('T')
    total = gc_count + at_count
    if total > 0:
        return (gc_count / total) * 100.0, gc_count, at_count
    return 0.0, 0, 0

def process_crop(crop_path, element_a, element_b):
    fasta_files = glob.glob(os.path.join(crop_path, "*.fasta"))
    if not fasta_files:
        return None
    
    total_spacer_g_c = 0
    total_spacer_a_t = 0
    genome_g_c = 0
    genome_a_t = 0
    
    for fasta_file in fasta_files:
        with open(fasta_file, 'r') as f:
            sequence = []
            for line in f:
                line = line.strip()
                if not line.startswith('>'):
                    sequence.append(line.upper())
            full_seq = "".join(sequence)
            
            _, g_c, a_t = calculate_gc(full_seq)
            genome_g_c += g_c
            genome_a_t += a_t
            
            promoters = find_promoters(full_seq)
            for prom in promoters:
                spacers = find_spacers(prom, element_a, element_b)
                for spacer in spacers:
                    _, g_c, a_t = calculate_gc(spacer)
                    total_spacer_g_c += g_c
                    total_spacer_a_t += a_t
                    
    total = total_spacer_g_c + total_spacer_a_t
    if total > 0:
        gc_pct = (total_spacer_g_c / total) * 100.0
    else:
        gc_pct = 0.0
        
    genome_total = genome_g_c + genome_a_t
    if genome_total > 0:
        genome_baseline = (genome_g_c / genome_total) * 100.0
    else:
        genome_baseline = 0.0
        
    return {
        "Spacer_GC_%": gc_pct,
        "Total_Spacers_Nucleotides": total,
        "Genome_Baseline_%": genome_baseline
    }

def main():
    print("=== Promoter Spacer GC Analysis ===")
    element_a = input("Enter first Cis-element: ").strip().upper()
    element_b = input("Enter second Cis-element (if same as first, enter it again): ").strip().upper()
    
    base_dir = "Data"
    if not os.path.exists(base_dir):
        print(f"Error: {base_dir} not found.")
        return
        
    crops = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    
    results = []
    labels = []
    values = []
    baselines = []
    
    for crop in crops:
        crop_path = os.path.join(base_dir, crop)
        print(f"Processing {crop}...")
        
        crop_data = process_crop(crop_path, element_a, element_b)
        if crop_data is None:
            continue
            
        results.append({
            "Crop": crop,
            "Spacer_GC_%": round(crop_data["Spacer_GC_%"], 2),
            "Total_Nucleotides": crop_data["Total_Spacers_Nucleotides"],
            "Genome_Baseline_%": round(crop_data["Genome_Baseline_%"], 2)
        })
        labels.append(crop)
        values.append(crop_data["Spacer_GC_%"])
        baselines.append(crop_data["Genome_Baseline_%"])
        
    # Tabular
    df = pd.DataFrame(results)
    csv_filename = f"Spacer_GC_Analysis_{element_a}_{element_b}.csv"
    df.to_csv(csv_filename, index=False)
    
    # Plotting
    plt.figure(figsize=(10, 6))
    
    # Different colors for each crop using a colormap
    colors = plt.cm.get_cmap('tab10', len(labels))
    bars = plt.bar(labels, values, color=[colors(i) for i in range(len(labels))])
    
    for bar, baseline in zip(bars, baselines):
        x = bar.get_x()
        w = bar.get_width()
        plt.hlines(baseline, x, x + w, colors='red', linestyles='dotted')
    
    plt.ylim(0, 50) # strictly forced y limit equivalent to the other graphs
    plt.xlabel('Crop')
    plt.ylabel('Spacer GC Content (%)')
    plt.title(f'Spacer GC Content for Cis-element Pair: {element_a} - {element_b}')
    plt.xticks(rotation=45, ha='right')
    
    plt.tight_layout()
    plot_filename = f"Spacer_GC_Plot_{element_a}_{element_b}.png"
    plt.savefig(plot_filename)
    plt.close()
    
    print(f"\nAnalysis complete.")
    print(f"Saved tabular data to: {csv_filename}")
    print(f"Saved bar plot to: {plot_filename}")

if __name__ == "__main__":
    main()
