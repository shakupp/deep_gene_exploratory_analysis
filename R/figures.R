library(ggplot2)
library(ggrepel)

# Volcano plot

volcano_plot = function(data_,
                        title_, # User title
                        subtitle_, # user subtitle
                        log2FC_, # For the threshold
                        padj_, # For the threshold
                        color_) { # For the expression status 
  ggplot(data = data_,
         aes(
           x = log2FC,
           y = -log10(padj) ,
           col = expressionStatus,
           label = GeneName,
         )) +
    geom_point() +
    theme_minimal() +
    geom_vline(xintercept = c(-log2FC_, log2FC_),
               col = "brown") +
    geom_hline(yintercept = -log10(padj_), col = "brown") +
    scale_colour_manual(values = color_) +
    geom_text_repel(data = data_[data_$expressionStatus != "NO",],
                    aes(label = GeneName),
                    max.overlaps = 10) +
    labs(
      title = title_,
      subtitle = subtitle_,
      y = "-log10(padj)",
      x = "log2(Fold Change)"
    )
}

# MA plot

MA_plot = function(data_,
                   title_, # User title
                   subtitle_, # User Subtitle
                   color_) { # For the expression status 
  ggplot(data = data_, aes(
    x = log2(baseMean), # BaseMean in log2 for better visualisation.
    y = log2FC,
    col = expressionStatus,
    label = GeneName
  )) +
    geom_point() +
    theme_minimal() +
    scale_colour_manual(values = color_) +
    geom_hline(yintercept = 0, col = "brown") +
    labs(
      title = title_,
      subtitle = subtitle_,
      y = "log2(Fold Change)",
      x = "log2(Base Mean)"
    )
}