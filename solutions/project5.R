library("phangorn")

data = read.dna("/Users/dmitrav/Library/Mobile Documents/com~apple~CloudDocs/ETHZ/Statistical models course/Lecture\ 5/ParisRT.txt")

# pairwise distances between sequences (Kimura model)
distance = dist.dna(data)

# create an initial tree topology for the alignment
initial.topology = nj(distance)

# plotting topology
plot(initial.topology)

# fitting model to the data
fit = pml(initial.topology, data = phyDat(data), model = "K80")

# log-likelihood value
fit$logLik

# optimise Q
fit.q = optim.pml(fit, optQ = TRUE)

# lower triangular part of the rate matrix
fit.q$Q

# optimise Q, edges, topology
fit.qet = optim.pml(fit, optQ = TRUE, optEdge = TRUE, optNni = TRUE)

# optimised log-likelihood
fit.qet$logLik

# bootstrap: resampling alignments
bs <- bootstrap.pml(fit.qet, bs=100, optNni=TRUE)

# plot shows: Mme S is more likely to have infected the patient Mme S 
treeBS <- plotBS(fit.qet$tree, bs, type = "phylogram")





