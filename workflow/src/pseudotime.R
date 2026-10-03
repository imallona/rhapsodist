## Pseudotime as the shortest-path distance from a root cell over a
## k-nearest-neighbour graph. Needs igraph.

## Index of the cell closest to the centroid of the cluster with the highest
## mean root score, such as the expression of a marker of the earliest state.
root_cell <- function(embedding, clusters, root_score) {
    cluster_means <- tapply(root_score, clusters, mean)
    in_root <- which(clusters == names(which.max(cluster_means)))
    centred <- sweep(embedding[in_root, , drop = FALSE], 2,
                     colMeans(embedding[in_root, , drop = FALSE]))
    in_root[which.min(rowSums(centred^2))]
}

## nn_index and nn_distance are cells by k matrices of neighbour indices and
## distances; a cell listed as its own neighbour is ignored. Returns distances
## from the root scaled to [0, 1], NA for cells not connected to it.
graph_pseudotime <- function(nn_index, nn_distance, root) {
    n_cells <- nrow(nn_index)
    edges <- data.frame(from = rep(seq_len(n_cells), ncol(nn_index)),
                        to = as.vector(nn_index),
                        weight = as.vector(nn_distance))
    edges <- edges[edges$from != edges$to, ]
    graph <- igraph::graph_from_data_frame(
        edges, directed = FALSE, vertices = data.frame(name = seq_len(n_cells)))
    from_root <- as.vector(igraph::distances(graph, v = root,
                                             weights = igraph::E(graph)$weight))
    from_root[is.infinite(from_root)] <- NA
    from_root / max(from_root, na.rm = TRUE)
}

## Spearman correlation of pseudotime for every pair of columns, each pair on
## the rows where both have a value.
pairwise_pseudotime_correlation <- function(pseudotime) {
    cor(pseudotime, method = "spearman", use = "pairwise.complete.obs")
}
