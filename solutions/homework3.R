#' Extract a specific emission probability from \data{emission}.
#'
#' @param state The state of the hidden variable.
#' @param observed The observed symbol.
#' @return The probability of `state` emitting `observed`.
#' @author Mathias Cardner
getEmission <- function(state, observed) {
  emission[state, observed]
}


#' Compute the forward quantity (message) at a given locus `n`.
#'
#' @param observed The observed symbol.
#' @param previous The forward quantity at locus `n-1`.
#' @return The forward quantity at the given locus.
#' @author Mathias Cardner
forward <- function(observed, previous) {
  # Compute the forward quantity using a matrix multiplication formulation.
  raw <- getEmission(c("+", "-"), observed) * transition %*% previous
  # Scale the result in order to prevent underflow.
  raw / sum(raw)
}

#' Compute the backward quantity (message) at a given locus `n`.
#'
#' @param observed The observed symbol.
#' @param upcoming The backward quantity at locus `n+1`.
#' @return The forward quantity at the given locus.
#' @author Mathias Cardner
backward <- function(observed, upcoming) {
  # Compute the backward quantity using a matrix multiplication formulation.
  raw <- transition %*% (getEmission(c("+", "-"), observed) * upcoming)
  # Scale the result in order to prevent underflow.
  raw / sum(raw)
}


# create a matrix for forward algorithm results
forward.values = matrix(rep(0, 2 * length(x)), nrow=2, ncol = length(x))
previous = c(0.5, 0.5)  # initial state

# compute other values with a forward algorithm
for (i in 1:length(x)){
  forward.values[,i] = forward(x[i], previous)
  previous = forward.values[,i]
}

# create a matrix for forward algorithm results
backward.values = matrix(rep(0, 2 * length(x)), nrow=2, ncol = length(x))
upcoming = c(1,1)  # initial state

# compute other values with a backward algorithm
for (j in length(x):1){
  backward.values[,j] = backward(x[j], upcoming)
  upcoming = backward.values[,j]
}

# multiply results
multiplied.values = forward.values * backward.values
# scale results
scaled.multiplied.values = apply(multiplied.values, 2, function(x){ x / sum(x) })

plot(Z == "+")
lines(Z == "+" & x == "C", col="lightblue")

