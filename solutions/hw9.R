
data = readRDS("../data/MVN_DAG.rds")

# plotting shows no correlation
attach(data)
plot(A, B, main="Scatterplot", xlab="A", ylab="B")

# # no correlation agrees with the grph structure, since A | B

# testing correlations
correlation.results = cor.test(data$A, data$B)

# test confirms there's no correlation
correlation.results$estimate
correlation.results$p.value


# regression A to C
ac.regression = lm(A ~ C, data = data)
ac.residuals = residuals(ac.regression)

# regression B to C
bc.regression = lm(B ~ C, data = data)
bc.residuals = residuals(bc.regression)

# plotting residuals
plot(ac.residuals, bc.residuals, main="Residuals", xlab="A ~ C", ylab="B ~ C")

# correlation test for residuals
cor.test(ac.residuals, bc.residuals)

# # there's negative correlation of residuals, which agrees wit hthe expectation,
# # since A and B are not independent given C


# PC algorithm
library(pcalg)

suffStat = list(C = cor(data), n = nrow(data))

pc.results.1 = pc(suffStat = suffStat, indepTest = gaussCItest, alpha = 1, labels = colnames(data), verbose=TRUE)
pc.results.2 = pc(suffStat = suffStat, indepTest = gaussCItest, alpha = 0.5, labels = colnames(data), verbose=TRUE)
pc.results.3 = pc(suffStat = suffStat, indepTest = gaussCItest, alpha = 0.2, labels = colnames(data), verbose=TRUE)
pc.results.4 = pc(suffStat = suffStat, indepTest = gaussCItest, alpha = 0.01, labels = colnames(data), verbose=TRUE)

par(mfrow=c(2,2))
plot(pc.results.1, main="Alpha = 1")
plot(pc.results.2, main="Alpha = 0.5")
plot(pc.results.3, main="Alpha = 0.25")
plot(pc.results.4, main="Alpha = 0.01")