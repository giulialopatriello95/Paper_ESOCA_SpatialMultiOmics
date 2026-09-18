############################################################
# Bulk RNA-seq analysis
# EAC response analysis
#
# Description:
# Differential expression analysis of bulk RNA-seq data,
# performed separately for the OSR and FREGAT cohorts.
# Results are subsequently combined by meta-analysis,
# followed by pathway enrichment analysis and PCA.
############################################################


############################
# 1. Load libraries
############################

library(org.Hs.eg.db)
library(DESeq2)
library(ggplot2)
library(edgeR)
library(SurfR)


############################
# 2. OSR cohort
############################

# Load gene-level expression matrix and sample metadata.
merged_df <- read.table(
  "aw_BulkRNAseq_ESOCA.xlsx",
  header = TRUE
)

sampledata <- read.table("sampledata.txt")


# Define the sequencing library type.
# Samples generated using Smart-seq are explicitly labelled.
sampledata$Library <- "mRNA"

sampledata$Library[
  sampledata$Sample %in% c(
    "RNAESOCA70T",
    "RNAESOCA77T",
    "RNAESOCA93T",
    "RNAESOCA97T",
    "RNAESOCA102T",
    "RNAFT9_bis"
  )
] <- "smartSeq"


# Retain only samples from the OSR cohort.
sampledata <- sampledata[sampledata$Cohort != "Fregat", ]


# Match the expression matrix to the samples present
# in the metadata.
merged_df <- merged_df[, rownames(sampledata)]


# Create the DESeq2 dataset.
# Differential expression is modelled according to treatment response.
dds <- DESeqDataSetFromMatrix(
  countData = merged_df,
  colData = sampledata,
  design = ~ Response
)

# Run the DESeq2 differential expression workflow.
dds <- DESeq(dds)

# Extract the normalized expression matrix.
mat <- assay(dds)


# Determine the smallest number of biological replicates
# among the response groups.
Nreplica_value.OSR <- min(table(sampledata$Response))


# Differential expression analysis using SurfR.
#
# CR = complete responder
# NR = non-responder
#
# Run is included in the design to account for technical variation.
df.OSR <- DGE(
  expression = merged_df,
  metadata = sampledata,
  Nreplica = Nreplica_value.OSR,
  design = "~ Run + Response",
  condition = "Response",
  alpha = 0.05,
  TEST = "CR",
  CTRL = "NR",
  output_tsv = FALSE
)


############################################################
# 3. FREGAT cohort
############################################################

# Load libraries required for the FREGAT analysis.
library(DESeq2)
library(ggplot2)


# Load expression matrix and metadata.
merged_df <- read.table(
  "aw_BulkRNAseq_ESOCA.xlsx",
  header = TRUE
)

sampledata <- read.table("sampledata.txt")


# Define sequencing library type.
sampledata$Library <- "mRNA"

sampledata$Library[
  sampledata$Sample %in% c(
    "RNAESOCA70T",
    "RNAESOCA77T",
    "RNAESOCA93T",
    "RNAESOCA97T",
    "RNAESOCA102T",
    "RNAFT9_bis"
  )
] <- "smartSeq"


# Retain only samples belonging to the FREGAT cohort.
sampledata <- sampledata[sampledata$Cohort == "Fregat", ]


# Match expression data to the selected samples.
merged_df <- merged_df[, rownames(sampledata)]


# Create and run the DESeq2 analysis.
dds <- DESeqDataSetFromMatrix(
  countData = merged_df,
  colData = sampledata,
  design = ~ Response
)

dds <- DESeq(dds)

# Extract normalized expression values.
mat <- assay(dds)


# Determine the minimum number of replicates
# among the response groups.
Nreplica_value.FREGAT <- min(table(sampledata$Response))


# Differential expression analysis using SurfR.
df.FREGAT <- DGE(
  expression = merged_df,
  metadata = sampledata,
  Nreplica = Nreplica_value.FREGAT,
  design = "~ Run + Response",
  condition = "Response",
  alpha = 0.05,
  TEST = "CR",
  CTRL = "NR",
  output_tsv = FALSE
)


############################################################
# 4. Meta-analysis of OSR and FREGAT
############################################################

# Combine differential expression results from the two cohorts
# using Fisher's method and the inverse-normal method.

L_fishercomb <- metaRNAseq(
  ind_deg = list(
    FREGAT = df.FREGAT,
    OSR = df.OSR
  ),
  test_statistic = "fishercomb",
  BHth = 0.05,
  adjpval.t = 0.05,
  nrep = c(
    Nreplica_value.FREGAT,
    Nreplica_value.OSR
  )
)


L_invnorm <- metaRNAseq(
  ind_deg = list(
    FREGAT = df.FREGAT,
    OSR = df.OSR
  ),
  test_statistic = "invnorm",
  BHth = 0.05,
  adjpval.t = 0.05,
  nrep = c(
    Nreplica_value.FREGAT,
    Nreplica_value.OSR
  )
)


# Combine the Fisher and inverse-normal meta-analysis results.
metacomb <- combine_fisher_invnorm(
  ind_deg = list(
    FREGAT = df.FREGAT,
    OSR = df.OSR
  ),
  invnorm = L_invnorm,
  fishercomb = L_fishercomb,
  adjpval = 0.05
)


############################################################
# 5. Gene set enrichment analysis
############################################################

library(msigdbr)

# Extract OSR log2 fold changes and use them as the
# ranked gene list for GSEA.
lfc_vector <- metacomb$OSR_log2FC
names(lfc_vector) <- metacomb$GeneID

# Rank genes from highest to lowest log2 fold change.
lfc_vector <- sort(lfc_vector, decreasing = TRUE)


# Retrieve MSigDB C2 gene sets for Homo sapiens.
msig_h <- msigdbr(
  species = "Homo sapiens",
  category = "C2"
) %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::rename(
    ont = gs_name,
    symbol = gene_symbol
)


# Set seed for reproducibility.
set.seed(2020)

library(ggplot2)


# Perform GSEA for the OSR cohort.
#
# pvalueCutoff = 1 retains the complete enrichment result,
# which is subsequently filtered for the pathways of interest.
gse_OSR <- clusterProfiler::GSEA(
  lfc_vector,
  pvalueCutoff = 1,
  scoreType = "std",
  TERM2GENE = msig_h,
  nPermSimple = 100000
)


# Repeat the GSEA using the FREGAT log2 fold changes.
lfc_vector <- metacomb$FREGAT_log2FC
names(lfc_vector) <- metacomb$GeneID

lfc_vector <- sort(lfc_vector, decreasing = TRUE)


gse_FREGAT <- clusterProfiler::GSEA(
  lfc_vector,
  pvalueCutoff = 1,
  scoreType = "std",
  TERM2GENE = msig_h,
  nPermSimple = 100000
)


############################################################
# 6. Compare pathway enrichment between cohorts
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)


# Extract normalized enrichment scores (NES) and FDR values
# from the OSR GSEA results.
df1 <- gse_OSR@result %>%
  dplyr::select(ID, NES, p.adjust) %>%
  dplyr::rename(
    OSR_NES = NES,
    OSR_FDR = p.adjust
  )


# Extract NES and FDR values from the FREGAT analysis.
df2 <- gse_FREGAT@result %>%
  dplyr::select(ID, NES, p.adjust) %>%
  dplyr::rename(
    FREGAT_NES = NES,
    FREGAT_FDR = p.adjust
  )


# Merge pathway-level results from both cohorts.
merged <- full_join(df1, df2, by = "ID")


# Define the pathways selected for visualization.
gene_sets <- c(
  "KEGG_ETHER_LIPID_METABOLISM",
  "KEGG_VEGF_SIGNALING_PATHWAY",
  "KEGG_CELL_CYCLE",
  "WP_TCELL_RECEPTOR_SIGNALING",
  "REACTOME_SIGNALING_BY_INTERLEUKINS",
  "REACTOME_NEUTROPHIL_DEGRANULATION",
  "REACTOME_TCR_SIGNALING",
  "REACTOME_GPCR_SIGNALING",
  "KEGG_CHEMOKINE_SIGNALING_PATHWAY",
  "SARRIO_EPITHELIAL_MESENCHYMAL_TRANSITION_UP",
  "PID_AURORA_A_PATHWAY",
  "REACTOME_GPCR_LIGAND_BINDING"
)


# Keep only the selected pathways.
merged <- merged %>%
  filter(ID %in% gene_sets)


# Convert the data to long format for plotting.
plot_df <- merged %>%
  pivot_longer(
    cols = c(OSR_NES, FREGAT_NES),
    names_to = "Condition",
    values_to = "NES"
  ) %>%
  mutate(
    FDR = ifelse(
      Condition == "OSR_NES",
      OSR_FDR,
      FREGAT_FDR
    ),
    Significance = ifelse(FDR < 0.05, "*", "")
  )


# Order pathways according to their OSR NES.
plot_df$ID <- factor(
  plot_df$ID,
  levels = plot_df %>%
    filter(Condition == "OSR_NES") %>%
    arrange(NES) %>%
    pull(ID)
)


# Remove database prefixes from pathway names.
plot_df$ID <- gsub("REACTOME_", "", plot_df$ID)
plot_df$ID <- gsub("KEGG_", "", plot_df$ID)
plot_df$ID <- gsub("SARRIO_", "", plot_df$ID)

plot_df$Condition <- gsub("_NES", "", plot_df$Condition)


# Classify pathways according to the response group
# associated with their enrichment.
plot_df$Enriched <- "CR enriched"

plot_df$Enriched[
  plot_df$ID %in% c(
    "PID_AURORA_A_PATHWAY",
    "EPITHELIAL_MESENCHYMAL_TRANSITION_UP",
    "ETHER_LIPID_METABOLISM",
    "VEGF_SIGNALING_PATHWAY",
    "CELL_CYCLE"
  )
] <- "NR enriched"


# Split pathways into CR- and NR-enriched groups.
CR <- plot_df[plot_df$Enriched == "CR enriched", ]
NR <- plot_df[plot_df$Enriched == "NR enriched", ]


# Define cohort order for plotting.
CR$Condition <- factor(
  CR$Condition,
  levels = c("OSR", "FREGAT")
)

NR$Condition <- factor(
  NR$Condition,
  levels = c("OSR", "FREGAT")
)


############################################################
# 7. Plot CR-enriched pathways
############################################################

p1 <- ggplot(
  CR,
  aes(
    x = Condition,
    y = ID,
    fill = NES
  )
) +
  geom_tile(color = "grey70") +
  geom_text(
    aes(label = Significance),
    color = "black",
    size = 5
  ) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0,
    name = "NES"
  ) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 14,
      family = "Arial",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      family = "Arial",
      color = "black"
    ),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "",
    x = "",
    y = "Pathway"
  ) +
  ggtitle("CR enriched\nsignatures")


############################################################
# 8. Plot NR-enriched pathways
############################################################

p2 <- ggplot(
  NR,
  aes(
    x = Condition,
    y = ID,
    fill = NES
  )
) +
  geom_tile(color = "grey70") +
  geom_text(
    aes(label = Significance),
    color = "black",
    size = 5
  ) +
  scale_fill_gradient2(
    low = "#4575b4",
    mid = "white",
    high = "#d73027",
    midpoint = 0,
    name = "NES"
  ) +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 14,
      family = "Arial",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 12,
      family = "Arial",
      color = "black"
    ),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    title = "",
    x = "",
    y = "Pathway"
  ) +
  ggtitle("NR enriched\nsignatures")


# Display the two pathway heatmaps together.
library(patchwork)

p1 + p2


############################################################
# 9. PCA analysis using DESeq2
############################################################

# Load libraries required for PCA and visualization.
library(org.Hs.eg.db)
library(DESeq2)
library(ggplot2)
library(edgeR)
library(SurfR)


# Load the bulk RNA-seq count matrix and metadata.
merged_df <- read.table(
  "aw_BulkRNAseq_ESOCA.xlsx",
  header = TRUE
)

sampledata <- read.table("sampledata.txt")


# Define sequencing library type.
sampledata$Library <- "mRNA"

sampledata$Library[
  sampledata$Sample %in% c(
    "RNAESOCA70T",
    "RNAESOCA77T",
    "RNAESOCA93T",
    "RNAESOCA97T",
    "RNAESOCA102T",
    "RNAFT9_bis"
  )
] <- "smartSeq"


# Match expression data to metadata.
merged_df <- merged_df[, rownames(sampledata)]


# Build the DESeq2 object.
# Batch is included in the design to account for
# potential technical effects.
dds <- DESeqDataSetFromMatrix(
  countData = merged_df,
  colData = sampledata,
  design = ~ Batch + Response
)

dds <- DESeq(dds)


# Extract differential expression results.
res <- results(dds)


# Transform the data using rlog for PCA visualization.
rld <- rlog(dds)


# Select genes with nominal p-value < 0.01.
genes <- res[
  res$pvalue < 0.01 &
    !is.na(res$pvalue < 0.01),
]

genes_use <- rownames(genes)


############################################################
# 10. PCA by treatment response
############################################################

# Generate PCA using the selected genes
# and colour samples according to response.
DESeq2::plotPCA(
  rld[genes_use, ],
  intgroup = c("Response")
) +
  theme(aspect.ratio = 1)


############################################################
# 11. PCA by cohort
############################################################

# Standardize the cohort name before plotting.
rld$Cohort <- gsub(
  "OSR_old",
  "OSR",
  rld$Cohort
)


# Generate PCA coloured by cohort.
p <- DESeq2::plotPCA(
  rld[rownames(genes), ],
  intgroup = c("Cohort")
) +
  scale_color_manual(
    values = c("#F54927", "#2A593A")
  ) +
  theme(aspect.ratio = 1) +
  geom_point(size = 5)


# Save the cohort PCA as a high-resolution JPEG.
ggsave(
  filename = "PCA_DESeq2_Cohort.jpg",
  plot = p,
  device = "jpg",
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)


############################################################
# 12. PCA by treatment response
############################################################

# Generate the response-based PCA.
p <- DESeq2::plotPCA(
  rld[rownames(genes_use), ],
  intgroup = c("Response")
) +
  theme(aspect.ratio = 1) +
  geom_point(size = 5)


# Save the response PCA.
ggsave(
  filename = "PCA_DESeq2_Response.jpg",
  plot = p,
  device = "jpg",
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)


############################################################
# End of analysis
############################################################
