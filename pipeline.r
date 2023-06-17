# Load the required packages
library(magrittr)
library(readr)
library(limma)
library(affy)
library(annotate)
library(pheatmap)
library(mouse4302.db)
library(affyPLM)
library(genefilter)

# Read the data
microarray_data <- ReadAffy()

# Read the targets file containing sample info
targets <- read_csv("targets.csv")
rownames(targets) <- targets$Name

sampleNames(microarray_data) <- sub(
  "\\.CEL\\.gz",
  "",
  sampleNames(microarray_data)
)

ordered_inds <- match(targets$Filename, sampleNames(microarray_data))
sampleNames(microarray_data) <- targets$Name[ordered_inds]

# Assign the phenotypic data
meta_data <- data.frame(labelDescription = c(
  "Sample Name Shorthand",
  "Filename Prefix", "Cell Line",
  "Experimental Treatment"
))

# Retain the rownames in targets[ordered_inds, ]
sorted_targets <- targets[ordered_inds, ]
rownames(sorted_targets) <- targets$Name[ordered_inds]

phenoData(microarray_data) <- new("AnnotatedDataFrame",
  data = sorted_targets,
  varMetadata = meta_data
)

# Save the targets file as .Rdata for Shiny app
save(targets, file = "targets.Rdata")

# Normalise with RMA -> outputs an ExpressionSet object, which is part of the
# Biobase OOP system; contains the expression matrix and phenotypic data
# Then filter to remove probes with low variance
normalised_data <- nsFilter(
  rma(microarray_data),
  remove.dupEntrez = FALSE,
  var.cutof = 0.5
)$eset

# Get the expression matrix
normalised_expr_matrix <- exprs(normalised_data)

# Create a lookup table mapping probe IDs to more useful gene IDs
gene_ids_lookup <- select(mouse4302.db,
  keys = rownames(normalised_data),
  columns = c("ENTREZID", "ENSEMBL", "SYMBOL"),
  keytype = "PROBEID"
)

# Find the index of each row of the expression set in gene_ids_lookup
lookup_inds <- match(rownames(normalised_data), gene_ids_lookup$PROBEID)

# Use the index to set the phenotypic data in the ExpressionSet
fData(normalised_data) <- gene_ids_lookup[lookup_inds, ]

# Find all rows that don’t have an EntrezID and remove then
valid_rows <- !is.na(fData(normalised_data)$ENTREZID)
normalised_data <- normalised_data[valid_rows, ]

# Build the design matrix
design_mat <- model.matrix(~ -1 + factor(c(1, 1, 1, 2, 2, 2)))
colnames(design_mat) <- c("Hypoxic", "Normoxic")

# This instructs Limma which comparisons to make
contrast_mat <- makeContrasts(
  Hypoxic - Normoxic,
  levels = design_mat
)

# Fit the model, make comparisons and get the results
summary_table <- lmFit(normalised_data, design_mat) %>%
  contrasts.fit(contrast_mat) %>%
  eBayes() %>%
  topTable(
    coef = 1,
    adjust = "fdr",
    number = nrow(normalised_data)
  )

summary_table$signif5pc <- summary_table$adj.P.Val < 0.05
summary_table$signif1pc <- summary_table$adj.P.Val < 0.01

# Save the limma results as .Rdata for Shiny app
save(summary_table, file = "summary_table.Rdata")

# Get the positions of probes from the limma results in the annotation dataframe
summary_inds <- match(summary_table$PROBEID, gene_ids_lookup$PROBEID)

# Append the entrez IDs to the limma results
summary_table$EntrezID <- gene_ids_lookup$ENTREZID[summary_inds]

pval_sorted_inds <- order(summary_table$adj.P.Val, decreasing = FALSE)
sorted_ids <- summary_table$PROBEID[pval_sorted_inds]

# Save the input data for Shiny app
sorted_exprs_matrix <- normalised_expr_matrix[sorted_ids, ]

save(sorted_exprs_matrix, file = "sorted_exprs_matrix.Rdata")

# Load the MSigDB hallmark gene sets for mouse
ms_hallmark <- readRDS("Mm.h.all.v7.1.entrez.rds")

# Convert to indexes
hallmark_lookup_inds <- ids2indices(
  ms_hallmark,
  fData(normalised_data)$ENTREZID
)

# Calculate P values of enrichment for the hallmark gene sets with camera
geneset_pvals <- camera(
  normalised_data,
  index = hallmark_lookup_inds,
  design = design_mat, contrast = contrast_mat[, 1]
)

# Significantly upregulated gene sets in the filtered data
sig_genesets_up <- geneset_pvals[geneset_pvals$PValue < 0.05 &
  geneset_pvals$Direction == "Up", ]

# Significantly downregulated gene sets in the filtered data
sig_genesets_down <- geneset_pvals[geneset_pvals$PValue < 0.05 &
  geneset_pvals$Direction == "Down", ]

# Save the enrichment results as .Rdata for Shiny app
save(sig_genesets_up, file = "sig_genesets_up.Rdata")
save(sig_genesets_down, file = "sig_genesets_down.Rdata")

# True/False vectors, used for colouring volcano plots
oxphin <- summary_table$EntrezID %in%
  ms_hallmark$HALLMARK_OXIDATIVE_PHOSPHORYLATION
hypoxin <- summary_table$EntrezID %in%
  ms_hallmark$HALLMARK_HYPOXIA
tnfin <- summary_table$EntrezID %in%
  ms_hallmark$HALLMARK_TNFA_SIGNALING_VIA_NFKB
fatin <- summary_table$EntrezID %in%
  ms_hallmark$HALLMARK_FATTY_ACID_METABOLISM
glycoxin <- summary_table$EntrezID %in%
  ms_hallmark$HALLMARK_GLYCOLYSIS

# Save the vectors for the shiny app
save(oxphin, hypoxin, tnfin, fatin, file = "gsea_for_volc_r.Rdata")
