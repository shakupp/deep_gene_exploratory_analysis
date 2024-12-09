# Function to format the dataframe
formating_dataframe = function(df) {
  names(df) = c("GeneName", "ID", "baseMean", "log2FC", "pval", "padj")
  df$GeneName = formating_GeneName(df$GeneName)
  df$ID = formating_ID(df$ID)
  df$baseMean = formating_baseMean(df$baseMean)
  df$log2FC = formating_log2FC(df$log2FC)
  df$pval = formating_pval(df$pval)
  df$padj = formating_pval(df$padj)

  return(df)
}

# I chose to separate the formatting functions to allow future modifications.

formating_GeneName = function(column_) {
  return(as.character(column_))
}

formating_ID = function(column_) {
  return(as.character(column_))
}

formating_baseMean = function(column_) {
  return(as.numeric(column_))
}

formating_log2FC = function(column_) {
  return(as.numeric(column_))
}

formating_pval = function(column_) {
  return(as.numeric(column_))
}

formating_padj = function(column_) {
  return(as.numeric(column_))
}


# Function that determines the status of the results table based on the sliders for logFC and padj.
set_differential_expression = function(data_,
                                       log2FC_,
                                       padj_) {
  data_$expressionStatus = "NO"
  data_$expressionStatus[data_$log2FC > log2FC_ &
                           data_$padj < padj_] = "UP"
  data_$expressionStatus[data_$log2FC < -log2FC_ &
                           data_$padj < padj_] = "DOWN"
  return(data = data_)
}
