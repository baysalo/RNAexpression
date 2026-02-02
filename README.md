# Host-Pathogen Interaction Analysis: Bacterial Antagonism of Fungal Pathogens

The analysis integrates **Fungal Transcriptomics (RNA-seq)** and **Bacterial Genomics** to propose a mechanism of biological control.

## 📂 Repository Structure

* `Master_Analysis_Script.R`: The complete R pipeline to reproduce all figures.
* `Diff_genes_foranalysis.csv`: The list of differentially expressed fungal genes (Downregulated).
* `Bacterial_Weapons_List.csv`: The output of the genomic mining script (identified virulence factors).
* `Figures/`: High-resolution TIFF outputs of the analysis.

## 📊 Analysis Pipeline

The provided R script (`Master_Analysis_Script.R`) performs the following 4 steps:

### 1. Transcriptome Overview
* **Volcano Plot:** Visualizes the global shift in fungal gene expression under bacterial stress.
* **Heatmap:** Shows the consistent downregulation of key virulence factors across biological replicates.

### 2. Bacterial Genomic Forensics
* **Weapon Mining:** Scans the bacterial genome (FASTA) for known antifungal machinery, including:
    * Chitinases & Glucanases (Cell wall degradation)
    * Siderophores (Iron theft)
    * Type VI Secretion Systems (Contact-dependent toxins)
    * Secondary Metabolites (Antibiotics)

### 3. Interaction Network Modeling
* **Gene-Specific Attack Map:** A hierarchical network graph linking specific bacterial genes (the "weapons") to the fungal pathways they target (the "victims").

## 🚀 How to Run the Code

1.  Clone this repository or download the files.
2.  Open `Master_Analysis_Script.R` in RStudio.
3.  Ensure the `.csv` data files are in your working directory.
4.  Run the script.
    * *Note:* Step 3 will prompt you to select your bacterial FASTA file (e.g., `.faa` or `.fsa_aa`).

## 🛠 Dependencies

The script requires the following R packages:
* `pheatmap`, `EnhancedVolcano` (Visualization)
* `igraph`, `ggraph` (Network Analysis)
* `Biostrings` (Genomic Sequence Handling)
* `dplyr`, `ggplot2` (Data Manipulation)

## 📄 Citation

If you use this code or data, please cite:
> **[Baysal]**, et al. (2025). "Bacillus velezensis EU07 suppresses Fusarium graminearum via transcriptomic reprogramming." *.
