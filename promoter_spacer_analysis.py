import os
import glob
from collections import Counter
import statistics
import pandas as pd

def find_promoters(sequence):
    promoters = []
    search_start = 0
    seq_len = len(sequence)
    
    while True:
        n = sequence.find("TATA", search_start)
        if n == -1:
            break
        
        if n > 1000:
            # Check for ATG between N+28 and N+35
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
    spacer_lengths = []
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
        # We DO NOT cap space_len <= 30 here as per the new requirement!
        spacer_lengths.append(space_len)
        
        pos = h + len(b)
        if pos <= l: 
            pos = l + 1
            
    return spacer_lengths

def process_crop(crop_path, element_a, element_b):
    all_spacer_lengths = []
    total_promoters = 0
    fasta_files = glob.glob(os.path.join(crop_path, "*.fasta")) # includes chr1.fasta, etc
    
    if not fasta_files:
        return None
        
    for fasta_file in fasta_files:
        with open(fasta_file, 'r') as f:
            sequence = []
            for line in f:
                line = line.strip()
                if not line.startswith('>'):
                    sequence.append(line.upper())
            
            full_seq = "".join(sequence)
            
            promoters = find_promoters(full_seq)
            total_promoters += len(promoters)
            
            for prom in promoters:
                spacers = find_spacers(prom, element_a, element_b)
                all_spacer_lengths.extend(spacers)
                
    return {
        "num_promoters": total_promoters,
        "spacers": all_spacer_lengths
    }

def main():
    print("=== Entire Promoter Spacer Analysis ===")
    element_a = input("Enter first Cis-element: ").strip().upper()
    element_b = input("Enter second Cis-element (if same as first, enter it again): ").strip().upper()
    
    base_dir = "Data"
    
    if not os.path.exists(base_dir):
        print(f"Error: {base_dir} not found.")
        return
        
    crops = [d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d))]
    
    print(f"\nFound {len(crops)} crop directories in {base_dir}")
    print(f"Analyzing spacers between {element_a} and {element_b}...\n")
    
    results = []
    
    for crop in crops:
        crop_path = os.path.join(base_dir, crop)
        print(f"Processing {crop}...")
        
        crop_data = process_crop(crop_path, element_a, element_b)
        if crop_data is None:
            print(f"  -> No .fasta files found in {crop}. Skipping.")
            continue
            
        spacers = crop_data["spacers"]
        num_promoters = crop_data["num_promoters"]
        
        if not spacers:
            results.append({
                "Crop": crop,
                "Promoters_Found": num_promoters,
                "Total_Spacers": 0,
                "Max_Spacer_Length": "N/A",
                "Avg_Spacer_Length": "N/A",
                "Mode_Spacer_Length": "N/A",
                "Median_Spacer_Length": "N/A"
            })
            continue
            
        max_len = max(spacers)
        avg_len = sum(spacers) / len(spacers)
        median_len = statistics.median(spacers)
        
        # Calculate Mode
        counts = Counter(spacers)
        mode_len = counts.most_common(1)[0][0]
        
        results.append({
            "Crop": crop,
            "Promoters_Found": num_promoters,
            "Total_Spacers": len(spacers),
            "Max_Spacer_Length": max_len,
            "Avg_Spacer_Length": round(avg_len, 2),
            "Mode_Spacer_Length": mode_len,
            "Median_Spacer_Length": median_len
        })
        
    # Display and Save Results
    df = pd.DataFrame(results)
    print("\n\n=== Analysis Results ===")
    print(df.to_string(index=False))
    
    output_filename = f"SpacerAnalysis_{element_a}_{element_b}.csv"
    df.to_csv(output_filename, index=False)
    print(f"\nResults saved to {output_filename}")

if __name__ == "__main__":
    main()
