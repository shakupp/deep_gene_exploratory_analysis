plot_barplot <- function(res, pval_thresh) {
  df <- as.data.frame(res@result)
  df <- df[df$p.adjust < pval_thresh, ]
  
  if (inherits(res, "gseaResult")) {
    ggplot() + theme_void() + annotate("text", x = 1, y = 1, label = "Barplot not available for GSEA")
  } else if (nrow(df) > 0) {
    barplot(res, showCategory = 10) + ggtitle("Barplot")
  } else {
    ggplot() + theme_void() + annotate("text", x = 1, y = 1, label = "No enriched term")
  }
}

plot_dotplot <- function(res, pval_thresh) {
  df <- as.data.frame(res@result)
  df <- df[df$p.adjust < pval_thresh, ]
  
  if (nrow(df) > 0) {
    dotplot(res, showCategory = 10) + ggtitle("Dotplot")
  } else {
    ggplot() + theme_void() + annotate("text", x = 1, y = 1, label = "No enriched term")
  }
}

plot_ridgeplot <- function(res, pval_thresh) {
  df <- as.data.frame(res@result)
  df <- df[df$p.adjust < pval_thresh, ]
  if (inherits(res, "enrichResult")) {
    plot.new(); text(0.5, 0.5, "Ridgeplot not available for ORA", cex = 1.2)
  } else if (inherits(res, "gseaResult") && nrow(df) > 0) {
    ridgeplot(res, showCategory = 10)
  } else {
    plot.new(); text(0.5, 0.5, "No enriched term", cex = 1.2)
  }
}

plot_gseaplot <- function(res, pval_thresh) {
  df <- as.data.frame(res@result)
  df <- df[df$p.adjust < pval_thresh, ]
  if (!inherits(res, "gseaResult")) {
    plot.new(); text(0.5, 0.5, "GSEA Plot is only available for GSEA method", cex = 1.2)
  } else if (nrow(df) == 0) {
    plot.new(); text(0.5, 0.5, "No enriched term", cex = 1.2)
  } else {
    gseaplot2(res, geneSetID = df$ID[1])
  }
}

plot_emapplot <- function(res, pval_thresh) {
  df <- as.data.frame(res@result)
  df <- df[df$p.adjust < pval_thresh, ]
  if (nrow(df) > 1) {
    emapplot(pairwise_termsim(res), showCategory = 10)
  } else {
    plot.new(); text(0.5, 0.5, "Not enough enriched terms for enrichment map", cex = 1.2)
  }
}

plot_cnetplot <- function(res, pval_thresh) {
  df <- as.data.frame(res@result)
  df <- df[df$p.adjust < pval_thresh, ]
  if (nrow(df) > 0) {
    cnetplot(res, showCategory = 10, foldChange = NULL)
  } else {
    plot.new(); text(0.5, 0.5, "No enriched term for cnetplot", cex = 1.2)
  }
}

plot_upsetplot <- function(res, pval_thresh) {
  df <- as.data.frame(res@result)
  df <- df[df$p.adjust < pval_thresh, ]
  if (nrow(df) > 1) {
    upsetplot(res)
  } else {
    plot.new(); text(0.5, 0.5, "Not enough enriched terms for upset plot", cex = 1.2)
  }
}
