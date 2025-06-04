# pathway.R

is_colnames_true <- function(cols) {
  all(c("GeneName", "ID", "baseMean", "log2FC", "pval", "padj") %in% cols)
}

read_valid_data <- function(file) {
  ext <- tools::file_ext(file$name)
  if (!ext %in% c("csv", "txt")) return(NULL)
  
  df_try1 <- try(read.csv(file$datapath, header = TRUE), silent = TRUE)
  if (inherits(df_try1, "try-error") || !is_colnames_true(colnames(df_try1))) {
    df_try2 <- try(read.csv(file$datapath, header = TRUE, sep = ";"), silent = TRUE)
    if (!inherits(df_try2, "try-error") && is_colnames_true(colnames(df_try2))) {
      return(df_try2)
    } else {
      return(NULL)
    }
  } else {
    return(df_try1)
  }
}

filter_input_data <- function(df, method, direction) {
  if (method == "ORA") {
    if (direction == "Over") {
      df <- df[df$log2FC > 0, ]
    } else if (direction == "Under") {
      df <- df[df$log2FC < 0, ]
    }
  }
  df
}

run_pathway_analysis <- function(df, method, database) {
  converted_ids <- bitr(df$ID, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  df <- merge(df, converted_ids, by.x = "ID", by.y = "ENSEMBL")
  if (nrow(df) == 0) return(NULL)
  
  df <- df %>% group_by(ENTREZID) %>% slice_max(order_by = abs(log2FC), n = 1) %>% ungroup()
  gene_list <- df$log2FC
  names(gene_list) <- df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  if (length(gene_list) == 0 || all(is.na(gene_list))) return(NULL)
  
  if (method == "ORA") {
    if (database == "KEGG") {
      enrichKEGG(gene = names(gene_list), organism = "mmu")
    } else {
      enrichPathway(gene = names(gene_list), organism = "mouse")
    }
  } else {
    if (database == "KEGG") {
      gseKEGG(geneList = gene_list, organism = "mmu")
    } else {
      gsePathway(geneList = gene_list, organism = "mouse")
    }
  }
}
