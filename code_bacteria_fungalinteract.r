cleancode.The Master Script: Bacterial-Fungal Interaction Analysis
# ==============================================================================
# MASTER SCRIPT: HOST-PATHOGEN INTERACTION ANALYSIS
# ==============================================================================

# --- STEP 0: Install and Load Required Libraries ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
required_packages <- c("pheatmap", "EnhancedVolcano", "igraph", "ggraph", 
                       "dplyr", "ggplot2", "stringr", "Biostrings", "readr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ==============================================================================
# STEP 1: LOAD FUNGAL DATA & VOLCANO PLOT
# ==============================================================================
print("--- Step 1: Analyzing Fungal Transcriptome ---")

# Load your data (Assuming it is already in your environment as 'Diff_genes_foranalysis')
# If not, uncomment the line below:
# Diff_genes_foranalysis <- read_csv("Diff_genes_foranalysis.csv")

# 1. Generate Volcano Plot
p_volcano <- EnhancedVolcano(Diff_genes_foranalysis,
    lab = Diff_genes_foranalysis$Genes,
    x = 'Fold', y = 'FDR',
    title = 'Volcano Plot: Fungal Pathogen vs Bacteria',
    pCutoff = 0.05, FCcutoff = 1.5,
    pointSize = 2.0, labSize = 3.0,
    col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
    legendPosition = 'right')

# Save Volcano Plot
ggsave("Figure1_Volcano_Plot.tiff", plot = p_volcano, width=8, height=8, dpi=300, bg="white")
print("Figure 1 saved.")

# ==============================================================================
# STEP 2: EXPRESSION HEATMAP
# ==============================================================================
print("--- Step 2: Generating Heatmap ---")

# Prepare Data: Filter Top 50 Genes
heatmap_data <- as.data.frame(Diff_genes_foranalysis) %>%
  arrange(FDR) %>%
  head(50)

rownames(heatmap_data) <- heatmap_data$Genes

# Select only expression columns (Robust selection)
heatmap_matrix <- heatmap_data %>% 
  dplyr::select(starts_with("C"), starts_with("T"))

# Annotations (Assuming 18 samples: 9 Control, 9 Treated based on your columns)
sample_info <- data.frame(Condition = factor(c(rep("Control", 9), rep("Treated", 9))))
rownames(sample_info) <- colnames(heatmap_matrix)

# Save Heatmap as TIFF
tiff("Figure2_Heatmap.tiff", width=8, height=10, units="in", res=300, compression="lzw")
pheatmap(heatmap_matrix, 
         cluster_rows = TRUE, cluster_cols = TRUE, scale = "row",
         annotation_col = sample_info,
         annotation_colors = list(Condition = c(Control="#2ecc71", Treated="#e74c3c")),
         show_rownames = TRUE, fontsize_row = 7,
         main = "Top 50 Downregulated Genes")
dev.off()
print("Figure 2 saved.")

# ==============================================================================
# STEP 3: BACTERIAL WEAPON MINING (Genomics)
# ==============================================================================
print("--- Step 3: Mining Bacterial Genome (Select File Now) ---")

# Select File via Pop-up
fasta_file <- file.choose()
print(paste("Analyzing file:", fasta_file))

proteins <- readAAStringSet(fasta_file)
gene_names <- names(proteins)

# Define Keywords
keywords <- list(
  "Cell_Wall_Degradation" = c("chitinase", "glucanase", "cellulase", "protease", "lipase", "chitin binding"),
  "Toxin_Production"      = c("polyketide", "non-ribosomal", "peptide synthetase", "toxin", "hemolysin", "cytolysin"),
  "Iron_Theft"            = c("siderophore", "iron uptake", "ferrichrome", "enterobactin"),
  "Secretion_Systems"     = c("type VI secretion", "type III secretion", "T6SS", "T3SS", "virulence"),
  "Antifungal_Antibiotics"= c("fungicin", "bacteriocin", "zwittermicin", "kanosamine", "antibiotic")
)

# Scan Genome
results <- data.frame(GeneID=gene_names, Description=gene_names, Category="Other", stringsAsFactors=FALSE)
for (cat in names(keywords)) {
  pattern <- paste(keywords[[cat]], collapse = "|")
  hits <- grepl(pattern, results$Description, ignore.case=TRUE)
  results$Category[hits] <- cat
}
weapons_found <- results %>% filter(Category != "Other")

# Save List and Plot
write.csv(weapons_found, "Bacterial_Weapons_List.csv", row.names=FALSE)

p_arsenal <- ggplot(weapons_found, aes(x=Category, fill=Category)) +
  geom_bar() + theme_bw() +
  labs(title="Bacterial Antifungal Arsenal", y="Gene Count", x="") +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_brewer(palette="Set1")

ggsave("Figure3_Bacterial_Arsenal.tiff", plot=p_arsenal, width=8, height=6, dpi=300, bg="white")
print("Figure 3 saved.")

# ==============================================================================
# STEP 4: GENE-SPECIFIC ATTACK MAP (Interaction Network)
# ==============================================================================
print("--- Step 4: Building Interaction Network ---")

# 1. Prepare Bacterial Nodes (Top 3 per category)
bacterial_nodes <- weapons_found %>%
  group_by(Category) %>% slice(1:3) %>% ungroup() %>%
  mutate(ShortName = str_trunc(Description, 20, "right")) %>%
  transmute(Source=ShortName, Target=Category, Type="BactGene")

# 2. Prepare Fungal Nodes (Top 10 Downregulated)
top_fungal_genes <- Diff_genes_foranalysis %>%
  arrange(FDR) %>% head(10) %>%
  transmute(Source="Transcriptome Response", Target=Genes, Type="FungGene")

# 3. Define Mechanisms
core_interactions <- data.frame(
  Source = c("Cell_Wall_Degradation", "Toxin_Production", "Iron_Theft", "Secretion_Systems", "Antifungal_Antibiotics"),
  Target = c("Cell Wall Integrity", "Mitochondrial Function", "Iron Homeostasis", "Immune Signaling (MAPK)", "Translation/Ribosomes"),
  Type = "Mechanism"
)

# 4. Link Mechanisms to Transcriptome
fungal_hub_links <- data.frame(Source=core_interactions$Target, Target="Transcriptome Response", Type="HubLink")

# 5. Combine Edges
all_edges <- rbind(
  bacterial_nodes %>% dplyr::select(Source, Target),
  core_interactions %>% dplyr::select(Source, Target),
  fungal_hub_links %>% dplyr::select(Source, Target),
  data.frame(Source=top_fungal_genes$Source, Target=top_fungal_genes$Target)
)

# 6. Plot
g <- graph_from_data_frame(all_edges, directed=TRUE)
V(g)$type <- case_when(
  V(g)$name %in% bacterial_nodes$Source ~ "Bacterial Gene",
  V(g)$name %in% bacterial_nodes$Target ~ "Weapon Category",
  V(g)$name %in% core_interactions$Target ~ "Fungal Target",
  V(g)$name == "Transcriptome Response" ~ "Fungal Target",
  TRUE ~ "Fungal Gene"
)

p_net <- ggraph(g, layout='kk') + 
  geom_edge_link(color="grey60", arrow=arrow(length=unit(2,'mm')), end_cap=circle(3,'mm')) +
  geom_node_point(aes(color=type, size=type)) +
  scale_color_manual(values=c("Bacterial Gene"="#e74c3c", "Weapon Category"="#c0392b", 
                              "Fungal Target"="#2980b9", "Fungal Gene"="#3498db")) +
  scale_size_manual(values=c("Bacterial Gene"=3, "Weapon Category"=6, "Fungal Target"=6, "Fungal Gene"=3)) +
  geom_node_text(aes(label=name), repel=TRUE, size=3, fontface="bold") +
  theme_void() +
  labs(title="Host-Pathogen Interaction Map")

ggsave("Figure4_Interaction_Map.tiff", plot=p_net, width=12, height=10, dpi=300, bg="white")
print("Figure 4 saved. All analysis complete!")
