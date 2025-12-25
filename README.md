**In order to run the code**


Put your data in Data Folder save it by name of crop and name each chromosome .fasta file as chr1.fasta, chr2.fasta, etc
type command : make
followed by: ./main

To delete generated .txt file: rm -f *.txt



**CIS ELEMENTS ANALYSIS FINAL REPORT**
---------------------------------------------
Rohit Pagar (2021B1A72139G)
Jaiaditya Agarwal (2021B1A72138G)


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
All results obtained can be accessed here: https://docs.google.com/spreadsheets/d/1U9pNqueA3ETWZrH6jkYG14FoqxvgRJuFnKYz0AzlT_8/edit?gid=0#gid=0
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
     4. Chi-Square and Fisher's Exact Tests:
Statistical tests were performed to assess whether the distribution of motifs between genome and promoter regions was significantly different:
Chi-Square Test: Evaluates differences in expected vs. observed frequencies of motif occurrences. The chi-square test compares observed and expected frequencies to determine if there is a significant association between two categorical variables.
χ² = Σ ((Oᵢ - Eᵢ)² / Eᵢ)
Oi​: Observed frequency in the iii-th category
Ei​: Expected frequency in the iii-th category
∑: Summation over all categories.

Fisher's Exact Test: Provides exact p-values, especially useful for datasets with small expected frequencies. Fisher's exact test is used for analyzing 2×2 contingency tables, particularly with small sample sizes. 
p = [(a+b)! × (c+d)! × (a+c)! × (b+d)!] / [a! × b! × c! × d! × n!]
n: Grand total of all observations (a+b+c+d)
a, b, c, d: Observed frequencies in each cell



Category 1
Category 2
Category 3
Group 1
a
b
a+b
Group 2
c
d
c+d
Col Total
a+c
b+d
n (Total)





Further Analysis

We compile data for enrichment ratio,Chi-square p-value, Fischer test p-value obtained for each cis-element - crop pair in  tabular form.
We try onto group crops into sets to continue our Further Analysis.

Fig-Tabular  compilation of data
Observations drawn from the table

Notably, the Pigeon Pea exhibits the highest enrichment ratios across all cis-elements (e.g., 2.52 for "ACGT-AAAG", 2.6 for "ACGT-ACGT", and 2.55 for "AAAG-AAAG") with statistically significant p-values (0 or near-zero), indicating a strong association of these elements with Pigeon Pea. Similarly, Chickpea demonstrates a relatively high enrichment for "ACGT-ACGT" (2.37) and "AAAG-AAAG" (1.78), supported by low p-values. These observations suggest a stronger binding affinity or functional significance of these cis-elements in these crops compared to others.
The higher enrichment ratios in Pigeon Pea might be attributed to unique genomic features, regulatory requirements, or evolutionary adaptations favoring these motifs. The statistical significance across tests reaffirms the reliability of the enrichment patterns. For crops like Sorghum and Quinoa, the lower ratios and relatively higher p-values imply weaker or less specific associations with these elements. This variability could stem from differences in regulatory networks, crop-specific environmental interactions, or sampling bias in motif detection.


To analyze the data comprehensively, we performed k-means clustering, plotted a dendrogram, and created a heat map to uncover patterns and relationships between crops and cis-elements. K-means clustering was used to group crops and cis-elements with similar enrichment ratios and statistical significance, allowing identification of clusters that share regulatory or genomic features. The dendrogram, generated from hierarchical clustering, provided a visual representation of the similarity and hierarchical relationships, highlighting groups of crops or cis-elements with common patterns. The heat map offered an intuitive way to visualize the entire dataset, making it easy to spot trends, such as high enrichment ratios and significant p-values in specific crops (e.g., Pigeon Pea) or cis-elements. These analyses are crucial for identifying crop-specific or shared regulatory motifs, understanding functional roles of cis-elements, and guiding targeted research, such as validating key motifs or exploring crop-specific regulatory networks. This multi-faceted approach helps prioritize biological questions and refine hypotheses for further experimental validation.

To analyze the data, we wrote the code using Python libraries such as pandas, numpy, seaborn, matplotlib, and clustering tools from scipy and sklearn. The first step involved normalizing the dataset using the StandardScaler to ensure that all variables contributed equally to the clustering process, as their scales could differ. we performed hierarchical clustering using the Ward method with the linkage() function to minimize variance within clusters and plotted a dendrogram to visually represent the relationships and distances between crops based on their enrichment ratios. For k-means clustering, we specified n_clusters=3 to identify distinct groupings of crops, using 10 initializations (n_init=10) to ensure stability and a fixed random_state=0 for reproducibility. The resulting clusters were assigned to the crops for further interpretation. Finally, we plotted a heat map with clustering using sns.clustermap() to visually correlate patterns across crops and cis-elements. This approach allowed me to systematically group crops and elements, aiding in identifying regulatory similarities, visualizing trends, and prioritizing biological features for further research.
Results Obtained from the code



Conclusion
The clustering results reveal three distinct groups among the crops based on enrichment ratios. Cluster 1 includes Sorghum, Fonio, Cassava, and Quinoa, indicating a close similarity in their cis-element enrichment profiles. Cluster 2 contains only Pigeon Pea, suggesting it has unique enrichment characteristics, likely driven by the significantly higher enrichment ratios and p-values seen in the initial table. Cluster 3 is represented solely by Chickpea, which also demonstrates distinct patterns compared to the other crops, possibly due to its relatively high enrichment in specific cis-elements. The hierarchical dendrogram supports these findings by showing the distances between clusters, emphasizing the distinctiveness of Pigeon Pea and Chickpea. The heatmap further visualizes these patterns, making it easier to interpret the relationships.
These insights are significant as they highlight potential crop-specific regulatory mechanisms or transcriptional responses. Moving forward, this analysis can direct further studies to investigate the biological roles of these cis-elements in specific clusters. For example, Pigeon Pea and Chickpea could be prioritized for functional validation of unique motifs or regulatory pathways, while Cluster 1 crops might be explored for shared regulatory strategies. Additionally, integrating these clusters with genomic or phenotypic data could provide deeper insights into the functional impact of these elements.


