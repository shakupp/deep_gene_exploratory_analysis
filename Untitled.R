# Load required libraries
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)

mouse_data <- read.csv("/Users/lmaouloud/Desktop/exemple_2025.csv", sep = ";") 

# 2. Adjust analysis parameters to be more lenient
test_params <- list(
  go_pval = 0.1,  # Increased from 0.05 to capture more terms
  qvalueCutoff = 0.3,  # Increased from 0.2
  inputGO = "BP",  # Biological Process
  species = "Mouse",
  go_analysisMethodChoice = c("Over Representation Analysis (ORA)"),
  go_oraChoice = c("Over expressed DEG", "Under expressed DEG")
)

# 3. Perform the GO analysis with better parameters
original_gene_list <- mouse_data$log2FC
names(original_gene_list) <- mouse_data$ID
gene_list <- sort(na.omit(original_gene_list), decreasing = TRUE)

# Select significant genes (using adjusted p-value)
sig_genes_df <- mouse_data[mouse_data$padj < test_params$go_pval, ]
genes <- sig_genes_df$log2FC
names(genes) <- sig_genes_df$ID
genes <- na.omit(genes)

# Split into up/down regulated
up_genes <- names(genes)[genes > 1]  # Only consider genes with |Log2FC| > 1
down_genes <- names(genes)[genes < -1]

# Get organism database
go_organism <- "org.Mm.eg.db"
universe <- names(gene_list)

### Test ORA with improved parameters
if ("Over Representation Analysis (ORA)" %in% test_params$go_analysisMethodChoice) {
  analyze_up <- "Over expressed DEG" %in% test_params$go_oraChoice
  analyze_down <- "Under expressed DEG" %in% test_params$go_oraChoice
  
  # ORA on upregulated genes (Log2FC > 1)
  if (analyze_up && length(up_genes) > 0) {
    result_up <- enrichGO(
      gene = up_genes,
      universe = universe,
      OrgDb = get(go_organism),
      keyType = "ENSEMBL",
      readable = TRUE,
      ont = test_params$inputGO,
      pvalueCutoff = test_params$go_pval,
      qvalueCutoff = test_params$qvalueCutoff,
      minGSSize = 3,  # Reduced from default 10
      maxGSSize = 500
    )
    print("ORA Upregulated Results:")
    if (!is.null(result_up) && nrow(result_up) > 0) {
      print(head(result_up))
      # Enhanced visualization
      print(dotplot(result_up, showCategory=15, title="Upregulated GO Terms") + 
              theme(axis.text.y = element_text(size=8)))
    } else {
      print("No significant upregulated ORA results found")
      print(paste("Number of upregulated genes:", length(up_genes)))
      print("Try further relaxing p-value cutoff or using more genes")
    }
  }
  
  # ORA on downregulated genes (Log2FC < -1)
  if (analyze_down && length(down_genes) > 0) {
    result_down <- enrichGO(
      gene = down_genes,
      universe = universe,
      OrgDb = get(go_organism),
      keyType = "ENSEMBL",
      readable = TRUE,
      ont = test_params$inputGO,
      pvalueCutoff = test_params$go_pval,
      qvalueCutoff = test_params$qvalueCutoff,
      minGSSize = 3,  # Reduced from default 10
      maxGSSize = 500
    )
    print("ORA Downregulated Results:")
    if (!is.null(result_down) && nrow(result_down) > 0) {
      print(head(result_down))
      print(dotplot(result_down, showCategory=15, title="Downregulated GO Terms") + 
              theme(axis.text.y = element_text(size=8)))
    } else {
      print("No significant downregulated ORA results found")
      print(paste("Number of downregulated genes:", length(down_genes)))
    }
  }
}