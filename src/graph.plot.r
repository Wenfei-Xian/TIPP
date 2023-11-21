args <- commandArgs(TRUE)

# Install and load the igraph package
# install.packages("igraph")
library(igraph)

# Read the file into a dataframe for edges
edges_df <- read.table(args[1], header = FALSE, stringsAsFactors = FALSE, row.names = NULL)
colnames(edges_df) <- c("from", "to")

# Create the graph
g <- graph_from_data_frame(edges_df, directed = FALSE)

# Read the cluster information
clusters <- readLines(args[2])
cluster_colors <- rainbow(length(clusters))

# Create a vector to store the color for each node
node_colors <- rep("green", vcount(g))
for (i in 1:length(clusters)) {
  nodes_in_cluster <- unlist(strsplit(clusters[i], "\t"))
  # Ensure that the nodes are present in the graph
  nodes_in_cluster <- nodes_in_cluster[nodes_in_cluster %in% V(g)$name]
  node_colors[V(g)[nodes_in_cluster]] <- cluster_colors[i]
}

# Plot the graph with the specified layout and save to a PDF
l <- layout_with_fr(g)
pdf(args[3], width = 8, height = 6)
plot(g, layout = l, vertex.size = 2, vertex.label.cex = 0.05, edge.width = 0.1, vertex.frame.width = 0.1, vertex.color = node_colors)
dev.off()
