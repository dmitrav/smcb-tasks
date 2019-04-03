
setwd("/Users/andreidm/ETH/courses/smcb/solutions")


### In this implementation the single covariance matrix computed from initial dataset is used.
### Only means are updated.

load("../data/sampling.rda")

markov.blanket = function(node.index, adjacency.matrix){
  
  node.children = which(adjacency.matrix[node.index, ] == 1)
  node.parents = which(adjacency.matrix[ ,node.index] == 1)
  node.coparents = which(which(adjacency.matrix[ ,node.children] == 1) != node.index)
  
  unique(c(node.children, node.parents, node.coparents))
}

# means and covariances
means = apply(gene.expression, 2, mean)
covariances = cov(gene.expression)

# create matrix to store results
particles.matrix = matrix(0, nrow = 1, ncol = 5)
colnames(particles.matrix) = colnames(gene.expression)

# initialise matrix
x0 = c(means[1:2], 0, means[4], 0)
particles.matrix[1,] = x0

# define constants
inversed.sigma.z1 = solve(covariances[-3,-3])
inversed.sigma.z2 = solve(covariances[-5,-5])

sigma.i.z1 = covariances[3,-3]
sigma.i.z2 = covariances[5,-5]

sigma.not_i.z1 = covariances[-3,3]
sigma.not_i.z2 = covariances[-5,5]

sigma.ii.z1 = covariances[3,3]
sigma.ii.z2 = covariances[5,5]

# sample
for (i in 1:11000){
  
  particle = particles.matrix[i,]
  
  # update means
  if (i == 1){
    means = particles.matrix[1,]
  } else {
    means = particles.matrix[i-1,]
    means[3] = mean(particles.matrix[,3])
    means[5] = mean(particles.matrix[,5])
  }
  
  # compute new parameters for YER124C
  beta0 = (means[3] - sigma.i.z1 %*% inversed.sigma.z1 %*% means[-3])
  beta = inversed.sigma.z1 %*% sigma.not_i.z1
  sigma.squared = sigma.ii.z1 - sigma.i.z1 %*% inversed.sigma.z1 %*% sigma.not_i.z1
  
  # draw new value for YER124C
  particle[3] = rnorm(1, beta0 + particle[-3] %*% beta, sigma.squared)
  
  # update means
  means[3] = mean(c(particles.matrix[,3], particle[3]))
  
  # compute new parameters for YNL327W
  beta0 = (means[5] - sigma.i.z2 %*% inversed.sigma.z2 %*% means[-5])
  beta = inversed.sigma.z2 %*% sigma.not_i.z2
  sigma.squared = sigma.ii.z2 - sigma.i.z2 %*% inversed.sigma.z2 %*% sigma.not_i.z2
  
  # draw new value for YNL327W
  particle[5] = rnorm(1, beta0 + particle[-5] %*% beta, sigma.squared)
  
  # add particle to the matrix
  particles.matrix = rbind(particles.matrix, particle)

}

# isolating genes of interest from burn-in and other genes
particles.matrix = particles.matrix[-(1:1000), c(3,5)]

# plotting samples
library(RColorBrewer)
rf <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
r <- rf(32)

library(hexbin)
h1 = hexbin(particles.matrix)
plot(h1, colramp=rf)

# means differ
sampled.means = apply(particles.matrix, 2, mean)
initial.means = apply(gene.expression, 2, mean)[c(3,5)]


