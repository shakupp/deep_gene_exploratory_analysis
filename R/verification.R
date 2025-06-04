# verification.R

# Verify if the file is in CSV format
is_csv_format = function(file_extension) {
  # Check if the file extension is not CSV
  if (!file_extension %in% c("csv")) {
    shinyalert(title = "Please upload a CSV file", type = "error")
    return(FALSE)  # Return FALSE if the file is not CSV
  }
  return(TRUE)  # Return TRUE if the file is CSV
}

# Verify if column names are correct
is_colnames_true = function(colnames_) {
  required_cols = c("GeneName", "ID", "baseMean", "log2FC", "pval", "padj")
  colnames_lower = tolower(colnames_)
  required_cols_lower = tolower(required_cols)
  
  # Check if all required columns are present
  if (!all(required_cols_lower %in% colnames_lower)) {
    shinyalert(title = "Error in column names",
               text = "The file columns must be: GeneName, ID, baseMean, log2FC, pval, padj.",
               type = "error")
    return(FALSE)  # Return FALSE if columns are incorrect
  }
  return(TRUE)  # Return TRUE if columns are correct
}

# Check if there are any NA values in the dataframe
is_NA = function(df) {
  # Check if there are any missing values (NA) in the data
  if (any(is.na(df))) {
    shinyalert(title = "Data Error: Missing Values",
               text = "Please remove rows with NA values or fix them.",
               type = "error")
    return(FALSE)  # Return FALSE if there are NAs
  }
  return(TRUE)  # Return TRUE if no NAs
}