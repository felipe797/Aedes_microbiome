### Function for removing rare features from an abundance table


remove_rare <- function(table, cutoff, min_reads) {
    # Filter columns based on minimum reads
    filtered <- table[, colSums(table) > min_reads]
    
    # Initialize vector to store column indices of features kept after filtering by cutoff
    colkept <- c()
    
    # Calculate the cutoff value based on the provided cutoff percentage
    cutoff <- ceiling(cutoff * nrow(filtered))
    
    # Iterate over each column in the filtered table
    for (i in 1:ncol(filtered)) {
        # Count the number of non-zero values in the column
        nonzero_col <- length(which(filtered[, i] > 0))
        
        # Check if the number of non-zero values exceeds the cutoff
        if (nonzero_col > cutoff) {
            # Keep the column index if it meets the criteria
            colkept <- c(colkept, i)
        }
    }
    
    # Return the filtered table with the selected columns
    return(filtered[drop = FALSE, , colkept])
}
