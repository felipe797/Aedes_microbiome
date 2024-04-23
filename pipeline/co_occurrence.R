### This function calculates correlation between microbial taxa abundances
### and computes key measures for a graph of co-occurrence analysis


co_occurrence_analysis <- function(abundance_file, cor_cutoff, p_cutoff) {
  # Load abundance data
  Abu <- read.delim(abundance_file, header = TRUE, row.names = 1)
  Abu <- as.matrix(Abu)
  
  # Filter OTUs by occurrence frequency
  table <- Abu
  table[table > 0] <- 1
  table.generalist <- Abu[rowSums(table) >= 25, ]
  Abu <- table.generalist
  
  # Perform co-occurrence network analysis
  pattern <- igraph::co_occurrence_network(Abu, cor_cutoff, p_cutoff)
  
  # Calculate network topological properties
  c <- igraph::cluster_walktrap(pattern$graph1)
  modularity_score <- igraph::modularity(c)
  md <- igraph::modularity(pattern$graph1, igraph::membership(c), weights = NULL)
  cc <- igraph::transitivity(pattern$graph1, vids = NULL, weights = NULL)
  spl <- igraph::average.path.length(pattern$graph1, directed = FALSE, unconnected = TRUE)
  gd <- igraph::graph.density(pattern$graph1, loops = FALSE)
  nd <- igraph::diameter(pattern$graph1, directed = FALSE, unconnected = TRUE, weights = NA)
  node_degree <- igraph::degree(pattern$graph1, v = igraph::V(pattern$graph1), mode = "all")
  ad <- mean(node_degree)
  
  global_topology <- data.frame(Modularity = modularity_score,
                                Modularity_Density = md,
                                Clustering_Coefficient = cc,
                                Average_Path_Length = spl,
                                Graph_Density = gd,
                                Network_Diameter = nd,
                                Average_Degree = ad)
  
  # Calculate node topological features
  betweenness_centrality <- igraph::betweenness(pattern$graph1, v = igraph::V(pattern$graph1),
                                                directed = FALSE, weights = NA,
                                                nobigint = TRUE, normalized = FALSE)
  closeness_centrality <- igraph::closeness(pattern$graph1, vids = igraph::V(pattern$graph1),
                                            weights = NA, normalized = FALSE)
  node_transitivity <- igraph::transitivity(pattern$graph1, type = c("local"), vids = NULL,
                                            weights = NA)
  
  node_topology <- data.frame(Node_Degree = node_degree,
                               Betweenness_Centrality = betweenness_centrality,
                               Closeness_Centrality = closeness_centrality,
                               Node_Transitivity = node_transitivity)
  
  # Create an abundance table for OTUs present in the network
  my_list <- row.names(pattern$matrix.cor3)
  logical <- row.names(Abu) %in% my_list
  tab_subset <- subset(Abu, logical)
  
  # Return data frames containing global topology, node topology, and network OTU abundance
  return(list(Global_Topology = global_topology,
              Node_Topology = node_topology,
              Network_OTU_Abundance = tab_subset))
}
