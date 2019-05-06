
load("../data/yeastStorey.rda")

# split the data
X = data[,-1]
y = data[, 1]

# set random seed
set.seed(123)

# get indexes for training set
training.set.size = floor(0.7 * nrow(X))
training.set.indexes = sample(seq_len(nrow(X)), size = training.set.size)

# get training set
X.train = X[training.set.indexes,]
X.test = X[-training.set.indexes,]

# get test set
y.train = y[training.set.indexes]
y.test = y[-training.set.indexes]

# elastic net
library(glmnet)

results.grid = matrix(0, nrow=3, ncol=11)
rownames(results.grid) = c("alpha", "lambda", "cvm")

results.grid[1,] = seq(0, 1, 0.1)

for (i in 1:length(alphas)){
  res = cv.glmnet(x=as.matrix(X.train), y=y.train, nfolds = 10, alpha = results.grid[1,i])
  results.grid[2,i] = res$lambda.min
  results.grid[3,i] = min(res$cvm)
}

results.grid

# # assuming you wanted exactly this...
# best.fit = glmnet(as.matrix(X.train), y.train, alpha = 1, lambda = 0.08842389)
best.fit = glmnet(as.matrix(X.train), y.train, alpha = 1)

y.predicted = predict.glmnet(best.fit, as.matrix(X.test), type="response")

# cmv error vs lambda
plot(results.grid[3,], log(results.grid[2,]), xlab = "Log(lambda)", ylab = "CV mean error", main = "CV results", type = "b")

# trace curve of coefficients as a function of log λ
plot(y=best.fit$df, x=log(best.fit$lambda), type = "s", xlab = "Log(lambda)", ylab = "Number of non-zero coefficients")

# roc
library("pROC")
y.predicted = y.predicted[,ncol(y.predicted)]
roc(response=y.test, predictor = y.predicted, plot=TRUE)

# get all the coefficients
coefs = coef(best.fit)
# show chosen ones
coefs@Dimnames[[1]][as.vector(coefs[,ncol(coefs)] != 0)]



