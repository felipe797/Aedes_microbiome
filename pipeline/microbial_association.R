### This script performs a Kruskal-Wallis test on the community data and 
### classifies the significant phylotypes by rankking its associations

### The script uses a physeq object containing a OTU table and a metadata data frame for the metagenomic samples 

microbe_significance_selection <- function(physeq, norm.meta = FALSE, select.variables = NULL, grouping_column, pvalue.threshold = 0.05) {
  
  # Normalize metadata if required
  if (norm.meta) {
    get_vars <- get.num.variables(physeq)  # Get numerical variables from phyloseq object
    norm_variables <- log(get_vars$num.variables / rowSums(get_vars$num.variables))  # Normalize numerical variables
    meta_table <- select.vars(norm_variables, get_vars$notnum.variables, select.variables)  # Select variables for metadata
    sample_data(physeq) <- meta_table  # Update metadata in phyloseq object
  } else if (!norm.meta) {  # If normalization is not required
    abund_table <- otu_table(physeq)  # Get abundance table from phyloseq object
    otu_table(physeq) <- log((abund_table + 1) / (rowSums(abund_table) + dim(abund_table)[2]))  # Log-transform abundance data
  }
  
  # Perform Kruskal-Wallis test on abundance data
  abund_table <- otu_table(physeq)
  meta_table <- data.frame(sample_data(physeq))
  meta_table$Groups <- meta_table[, grouping_column]
  
  kruskal.wallis.table <- data.frame()
  data <- as.data.frame(abund_table)
  for (i in 1:dim(data)[2]) {
    ks.test <- kruskal.test(data[, i], g = meta_table$Groups)  # Perform Kruskal-Wallis test for each feature
    kruskal.wallis.table <- rbind(kruskal.wallis.table, data.frame(id = names(data)[i], p.value = ks.test$p.value))  # Store results
  }
  
  # Adjust p-values and select significant features
  kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]
  kruskal.wallis.table$FWER <- pbinom(q = 0, p = kruskal.wallis.table$p.value, size = dim(kruskal.wallis.table)[1], lower.tail = FALSE)
  kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value, decreasing = FALSE), ]
  kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
  kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor
  rownames(kruskal.wallis.table) <- kruskal.wallis.table$id
  
  last.significant.element <- max(which(kruskal.wallis.table$q.value <= pvalue.threshold))
  selected <- 1:last.significant.element
  sig_res <- kruskal.wallis.table$id[selected]
  
  # Prepare data for random forest classifier
  subset_data <- data.frame(data[, as.character(kruskal.wallis.table[rownames(kruskal.wallis.table), "id"])])
  kruskal.wallis.table$id <- colnames(subset_data)
  subset_data <- subset_data[, sig_res]
  
  # Perform random forest classifier
  IDs_map <- data.frame(row.names = colnames(subset_data), "taxa" = colnames(subset_data))
  val <- randomForest::randomForest(meta_table$Groups ~ ., data = subset_data, importance = TRUE, proximity = TRUE, ntree = 1500, keep.forest = FALSE)
  imp <- randomForest::importance(val)
  df_accuracy <- data.frame(row.names = NULL, Sample = rownames(imp), Value = abs(as.numeric(imp[, "MeanDecreaseAccuracy"])), Index = rep("Mean Decrease Accuracy", dim(imp)[1]))
  
  # Rearrange feature importance data for plotting
  df_accuracy$Sample <- IDs_map[as.character(df_accuracy$Sample), "taxa"]
  df_accuracy_order <- as.character(IDs_map[rownames(imp), "taxa"][order(abs(as.numeric(imp[, "MeanDecreaseAccuracy"])), decreasing = TRUE)])
  df_accuracy$Sample <- factor(as.character(df_accuracy$Sample), levels = df_accuracy_order)
  df_accuracy$rank <- base::rank(df_accuracy$Value, ties.method = "min")
  df_accuracy$rank <- max(df_accuracy$rank) - df_accuracy$rank + 1
  
  # Prepare data for plotting
  df <- NULL
  for (i in df_accuracy$Sample) {
    rank <- (subset(df_accuracy, df_accuracy$Sample == i))$rank
    tmp <- data.frame(subset_data[, i], meta_table$Groups, rep(rank), rep(paste(i, " p.adj = ", sprintf("%.10g", kruskal.wallis.table[kruskal.wallis.table$id == i, "q.value"]), sep = ""), dim(data)[1]))
    colnames(tmp) <- c("Value", "Groups", "Rank", "Taxa")
    if (is.null(df)) {
      df <- tmp
    } else {
      df <- rbind(df, tmp)
    }
    df <- na.omit(df)
  }
  
  # Output results as a list
  out <- list("SignfeaturesTable" = kruskal.wallis.table, "plotdata" = df, "importance" = df_accuracy)
  return(out)
}
