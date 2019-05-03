
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
?glmnet

