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


# NOT SURE AT ALL

# create a list for forward algorithm results and add the first value for initial state
forward.values = list(data.matrix(c(0.5,0.5)) %*% getEmission(c("+", "-"), x[1]) * transition)

# compute other values with a forward algorithm
for (i in 2:length(x)){
    forward.values[[i]] = forward(x[i], forward.values[[i-1]])
}

# create a list for forward algorithm results and add the first value for initial state
backward.values = list( transition %*% (data.matrix(c(1,1)) %*% getEmission(c("+", "-"), x[300])) )

# compute other values with a backward algorithm
for (i in 1:(length(x)-1) ){
  backward.values[[i+1]] = backward(x[length(x)-i], backward.values[[i]])
}

forward.values[[1]] %*% backward.values[[1]]

transition %*% (getEmission(c("+", "-"), x[299]) * upcoming)




