
setwd("/Users/andreidm/ETH/courses/smcb/solutions")

load("../data/sampling.rda")

markov.blanket = function(node.index, adjacency.matrix){
  
  node.children = which(adjacency.matrix[node.index, ] == 1)
  node.parents = which(adjacency.matrix[ ,node.index] == 1)
  node.coparents = which(adjacency.matrix[ ,node.children] == 1)
  
  unique(c(node.children, node.parents, node.coparents))
}

# means and covariances
means = apply(gene.expression, 2, mean)
covariances = cov(gene.expression)

# sampling
for (i in 1:11000){
  
}