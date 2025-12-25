**In order to run the code**


Put your data in Data Folder save it by name of crop and name each chromosome .fasta file as chr1.fasta, chr2.fasta, etc
type command : make
followed by: ./main

To delete generated .txt file: rm -f *.txt



**CIS ELEMENTS ANALYSIS FINAL REPORT**
---------------------------------------------

The primary objective of this study was to identify distinct patterns in spacer sequences between specific cis-element pairs and investigate their correlation with drought-related traits in various crops. To achieve this, we developed a comprehensive computational approach that systematically analyzes crop genomes on a chromosome-by-chromosome basis. This methodology identifies promoter regions, extracts spacer sequences ranging from 0 to 20 base pairs between targeted cis-element pairs, and generates chromosome-specific output in .txt format.
How the Code Works:
Input Parameters: The code takes three inputs: the cis-element pair to analyze, the crop genome name, and the number of chromosomes in the genome.
Spacer Sequence Extraction:
The program scans the entire genome to locate occurrences of the first cis-element in the pair. For each identified position i, it searches the subsequent sequence (i+1 to i+20) for the second cis-element in the pair.
If the second cis-element is found within this range, the spacer sequence (the sequence between the two elements) is added to a spacer vector.
If the second cis-element is not found, the algorithm moves to the next occurrence of the first cis-element and repeats the process until the entire genome is analyzed.
Promoter Sequence Identification:
In the same chromosome, the program searches for promoter sequences by identifying the presence of a TATA box in the genome.
For every TATA box found, it searches for an ATG sequence located 28-35 base pairs downstream.
If this condition is met, a region spanning 1000 base pairs upstream from the ATG sequence is extracted as a promoter sequence and added to the promoter vector.
The extracted promoter sequences are then subjected to the same spacer sequence extraction process, ensuring spacer sequences within promoter regions are thoroughly analyzed.
Chromosome-Wise Analysis:
This process is repeated for all chromosomes (i = 1 to n, where n is the total number of chromosomes).
For each chromosome, the spacer sequence data is saved into .txt files, organized chromosome by chromosome.
Data Compilation and Visualization:
Once the analysis is complete, data from all chromosomes is compiled into a consolidated spreadsheet.
A graph of Frequency vs. Spacer Length is plotted to visually identify any anomalies or distinct patterns that may emerge across different crops or cis-element pairs.
Crops and Cis-Elements Analyzed:
Monocots: Sorghum, Fonio
Dicots: Pigeon Pea, Quinoa, Cassava, Chickpea
Cis-Element Pairs: ACGT-ACGT, ACGT-AAAG, AAAG-AAAG
All results obtained can be accessed here: //link
Motif Occurrence Analysis
Following the identification and extraction of spacer sequences between specific cis-element pairs, our next step was to analyze the motif occurrences within the entire genome and promoter regions of the selected crops. This analysis aimed to evaluate the density, proportion, and enrichment of cis-elements, allowing us to explore their role in drought-related characteristics.

Methodology:
Motif Occurrence Analysis:
For each cis-element pair (ACGT-ACGT, ACGT-AAAG, AAAG-AAAG), we counted the total occurrences across the entire genome and within promoter regions for all selected crops.
The total motif occurrences for each genome and its promoters were recorded, enabling crop-to-crop comparisons.

Density Calculations:
Genome Density: For each motif, the density was calculated as the total occurrences divided by the genome size (in base pairs). 
Genome Density=Total Motif Occurrences in Genome/Genome Size ​
Promoter Density: Similarly, the density within promoters was computed as the total motif occurrences in promoters divided by the promoter size (in base pairs). Promoter Density=Total Motif Occurrences in Promoters/Promoter Size

Proportion of Motifs in Promoters: To assess the relative abundance of motifs in promoter regions compared to the genome, the proportion was calculated as: Proportion in Promoters=Total Motif Occurrences in Promoters/Total Motif Occurrences in Genome×100
Enrichment Ratio:
The enrichment ratio was determined to quantify the degree of overrepresentation (or underrepresentation) of motifs in promoter regions compared to the genome. Enrichment Ratio=Promoter Density/Genome Density











