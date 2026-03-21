##This program creates chromosome wise .fasta files from whole genome.fa/fasta files

from Bio import SeqIO
import re
import os

input_fasta = "./Data/Quinoa/quinoa_genome.fasta" ##(ALTER THIS)
output_dir = "./Data/Quinoa"  ##(ALTER THIS)
os.makedirs(output_dir, exist_ok=True)

for record in SeqIO.parse(input_fasta, "fasta"):
    
    # Combine id + description for maximum robustness
    text = record.id + " " + record.description
    
    # Try to extract chromosome number + A/B
    match = re.search(r'(?:chr|chromosome)[^\d]*(\d+)\s*([AB])', text, re.IGNORECASE)
    
    if match:
        chrom_num = int(match.group(1))     # e.g. 1
        subgenome = match.group(2).upper()  # A or B
        
        # Map to 1–18 numbering
        if subgenome == "A":
            new_index = 2 * chrom_num - 1
        else:
            new_index = 2 * chrom_num
        
        output_file = f"{output_dir}/chr{new_index}.fasta"
        
        with open(output_file, "w") as f:
            SeqIO.write(record, f, "fasta")
    
    else:
        print(f"⚠️ Skipping: {record.description}")

print("Done!")