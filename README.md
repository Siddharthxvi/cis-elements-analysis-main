# Cis-Regulatory Element (CRE) Architecture & Spacer Analysis Pipeline

A comprehensive, high-throughput computational pipeline for identifying, extracting, and analyzing the spatial and compositional architecture of *cis*-regulatory elements (CREs) across entire plant macro-genomes.

This repository focuses on evaluating the **geometric constraints (spacer distributions)** and **evolutionary compositional stability (GC content)** of co-occurring transcription factor binding sites within defined promoter micro-environments.

---

## 🧬 Overview

The spatial organization of *cis*-regulatory elements is fundamentally tied to transcriptional regulation. While single motifs are frequently studied, the syntax governing how two or more motifs interact—specifically the distance (spacer) and compositional environment between them—is a critical frontier in synthetic biology and evolutionary genomics.

This pipeline provides a robust C++ core for rapid, chromosome-scale sequence parsing combined with a suite of Python-based statistical modules to evaluate:
1. **Promoter Motif Co-occurrence:** Rapid extraction of 1000 bp upstream promoters strictly utilizing TATA box and downstream ATG anchoring.
2. **Spacer Distribution:** Unbiased calculation of spatial distribution spread (SD, IQR) without artificial boundary caps, explicitly evaluating structural fusion vs. stochastic distancing.
3. **Composite Element Identification:** Automated detection and filtering of 0-bp contiguous or overlapping motif boundaries (fused elements).
4. **Compositional Stability:** Parametric validation (ANOVA) and Coefficient of Variation (CV) tests measuring the stability of spacer GC content globally against raw genome-wide baselines.

> **Note on Sequence Logos:** The downstream generation of visual motif matrices and sequence logos for these extracted alignments is handled in our complementary repository. *(See: [Sequence Logo Analysis Repo](#related-repositories))*

---

## 🚀 Key Features and Modules

### 1. C++ Core Engine
The high-performance core iterates through massive FASTA arrays directory by directory.
*   **`Main.cpp` / `FindPromoter.cpp` / `FindSpacers.cpp`:** Reconstructs full chromosome sequences, anchors promoter boundaries using rigid `TATA -> ATG` positional rules, and systematically isolates nucleotide sequences occurring exclusively between Motif A and Motif B.

### 2. Spacer Distribution Statistics (`spacer_distribution_analysis.py`)
Calculates the degree of spatial variability between paired motifs.
*   Identifies robust central tendencies (Mean, Median, Mode).
*   Separates **composite elements** ($\le$ 0 bp spacer lengths) from natively spaced elements to present bias-free descriptive stats.
*   Outputs highly structured CSV arrays with IQR and SD metrics demonstrating physical constraint boundaries.

### 3. GC Content Profiling (`gc_content_analysis.py` & `spacer_gc_analysis.py`)
Tests for local evolutionary sequence constraints versus broad genomic mutational drift.
*   Computes exact genome-wide `%GC` to serve as a strict macro-baseline for every individual crop.
*   Overlays the local `%GC` composition of isolated *cis*-element spacers onto visual bar charts using the genome-baseline as a visual threshold.
*   Includes rigidly bounded (0-50%) multi-crop visualization graphics.

### 4. Chromosomal Stability Validation (`chromosomal_stability_stats.py`)
Mathematically proves spatial/chromosomal homogeneity.
*   **One-Way ANOVA:** Validates that variance in spacer GC content between different chromosomes within the same crop is statistically insignificant.
*   **Coefficient of Variation (CV):** Proves mathematical consistency (typically $< 2\%$) in nucleotide composition regardless of spatial chromosomal geography.

---

## 🛠 Prerequisites & Installation

### Environment Requirements
*   **C++ Compiler:** `g++` (Standard support for C++11 or later)
*   **Python:** `3.8+`
*   **Python Libraries:** `numpy`, `pandas`, `matplotlib`, `scipy`

### Setup
Clone the repository and ensure your biological sequence data is placed correctly:
```bash
git clone https://github.com/your-username/Cis-element-analysis-main.git
cd Cis-element-analysis-main

# Install required statistical packages
pip3 install numpy pandas matplotlib scipy
```

Ensure your genome files are organized natively inside the `Data/` folder by crop.
```text
Data/
├── Arabidopsis/
│   ├── chr1.fasta
│   └── chr2.fasta
├── Sorghum/
│   └── chr1.fasta
└── ...
```

---

## ⚙️ Usage Workflow

### Step 1: Extract Core Alignments (C++ Pipeline)
Compile the sequence parsing engine and run the binary to generate isolated promoter-FASTA files for specific pairs (e.g., `AAAG` and `ACGT`).
```bash
make clean && make
./Main
```
*(You will be prompted to enter your CREs, crop directory name, and starting chromosome).*

### Step 2: Distribution Analytics
To generate the raw spatial arrays and composite counts for all available crops simultaneously:
```bash
python3 spacer_distribution_analysis.py
```
*Outputs:* `Spacer_Distribution_Statistics.csv`

### Step 3: GC Composition Plots
To extract spacer sequences and plot them against the global organism baseline:
```bash
python3 spacer_gc_analysis.py
```
*(Enter the desired motif pair when prompted. The tool outputs a standardized `.png` graphic and CSV).*

### Step 4: Validate Chromosomal Stability
To calculate the ANOVA F-statistics, p-values, and CV mapping for your combinations:
```bash
python3 chromosomal_stability_stats.py
```
*Outputs:* `Chromosomal_Stability_Analysis.csv`

---

## 📊 Pipeline Parameters

The core configurable parameters that dictate boundary conditions within this pipeline:

| Parameter | Default Value | Description |
| :--- | :--- | :--- |
| **Promoter Window Size** | `1000 bp` | The upstream genomic region defined as the promoter (anchored strictly 28-35 nt upstream of an `ATG` downstream from a `TATA` box). |
| **Spacer Length Range** | Unbounded | The allowable spatial distance between two motif occurrences. Lengths of `0 bp` explicitly denote composite/adjacent motifs. |
| **Conservation Threshold** | `0.50` (50%) | Minimum frequency required for a nucleotide to be classified as positionally conserved inside generated consensus strings. |

---

## 🔗 Related Repositories
The downstream statistical mapping, thresholding, and visual generation of nucleotide **Sequence Logos** based on the FASTA alignments generated by this C++ core are housed separately. 

Please see the **Sequence Logo Analytics Repository** for generating quantitative PWM (Position Weight Matrix) visuals and thresholded sequence logos.

---

## 📝 Citation & Contact
If you utilize this pipeline or its statistical logic in your research or publication, please cite the associated manuscript (citation pending).

For issues, please open a ticket on GitHub or submit a Pull Request.
