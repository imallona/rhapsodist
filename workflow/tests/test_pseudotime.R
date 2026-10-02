#!/usr/bin/env Rscript
## Self-contained tests for workflow/src/pseudotime.R.
## Run with:
##   cd workflow && Rscript tests/test_pseudotime.R

script_dir <- local({
    ca <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", ca, value = TRUE)
    if (length(file_arg) > 0) {
        dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
        "."
    }
})
source(file.path(script_dir, "..", "src", "pseudotime.R"))

check <- function(label, condition) {
    if (!isTRUE(condition)) {
        stop(sprintf("FAIL: %s", label), call. = FALSE)
    }
    cat(sprintf("ok  - %s\n", label))
}

## k nearest neighbours by brute force, each cell first as its own neighbour
nearest_neighbours <- function(embedding, k) {
    d <- as.matrix(dist(embedding))
    index <- t(apply(d, 1, function(row) order(row)[seq_len(k + 1)]))
    distance <- t(sapply(seq_len(nrow(d)), function(i) d[i, index[i, ]]))
    list(index = index, distance = distance)
}

## 50 cells along a line, in shuffled order
set.seed(1)
position <- sample(seq(0, 10, length.out = 50))
line <- cbind(position, 0)
nn <- nearest_neighbours(line, k = 3)
root <- which.min(position)
pseudotime <- graph_pseudotime(nn$index, nn$distance, root)

check("the root has pseudotime 0", pseudotime[root] == 0)
check("the farthest cell has pseudotime 1", pseudotime[which.max(position)] == 1)
check("pseudotime follows the position on a line",
      isTRUE(all.equal(pseudotime, position / 10)))

## a second group far away, with neighbours only inside its own group
two_groups <- rbind(line, cbind(seq(100, 101, length.out = 5), 0))
nn_groups <- nearest_neighbours(two_groups, k = 3)
pseudotime_groups <- graph_pseudotime(nn_groups$index, nn_groups$distance, root)
check("cells not connected to the root get NA", all(is.na(pseudotime_groups[51:55])))
check("connected cells keep their pseudotime",
      isTRUE(all.equal(pseudotime_groups[1:50], pseudotime)))

clusters <- ifelse(position < 2, "early", ifelse(position < 7, "middle", "late"))
root_score <- 10 - position
picked <- root_cell(line, clusters, root_score)
check("the root is taken from the cluster with the highest score", clusters[picked] == "early")
early <- which(clusters == "early")
check("the root is the cell nearest the centroid of that cluster",
      picked == early[which.min(abs(position[early] - mean(position[early])))])

shifted <- pseudotime^2
shifted[1:5] <- NA
correlation <- pairwise_pseudotime_correlation(cbind(a = pseudotime, b = shifted, c = -pseudotime))
check("a monotone change keeps correlation 1 on shared cells",
      isTRUE(all.equal(correlation["a", "b"], 1)))
check("a reversed order gives correlation -1", isTRUE(all.equal(correlation["a", "c"], -1)))

cat("\nAll pseudotime tests passed.\n")
