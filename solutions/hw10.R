
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

for (i in 1:length(results.grid[1,])){
  res = cv.glmnet(x=as.matrix(X.train), y=y.train, nfolds = 10, family = "binomial", alpha = results.grid[1,i])
  results.grid[2,i] = res$lambda.min
  results.grid[3,i] = min(res$cvm)
}

results.grid

best.alpha = results.grid[1,results.grid[3,] == min(results.grid[3,])]
best.lambda = results.grid[2,results.grid[3,] == min(results.grid[3,])]


# ok, I'm getting the fit
best.fit = glmnet(as.matrix(X.train), y.train, family = "binomial", alpha = best.alpha, lambda = best.lambda)

# predict with "the best" model
y.predicted = predict.glmnet(best.fit, as.matrix(X.test), type="response")

# plot error vs log lambda
plot(cv.glmnet(as.matrix(X.train), y.train, alpha = best.alpha, family = "binomial"))

# trace curve of coefficients as a function of log lambda
plot(glmnet(as.matrix(X.train), y.train, alpha = best.alpha, family = "binomial"), xvar = "lambda")

# roc
library("pROC")
roc(response=y.test, predictor = y.predicted, plot=TRUE)

# get all the coefficients
coefs = coef(best.fit)

# show chosen ones
coefs@Dimnames[[1]][as.vector(coefs[,ncol(coefs)] != 0)][-1]



