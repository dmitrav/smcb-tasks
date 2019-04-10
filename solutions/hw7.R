
p = function(x1=NULL, x2=NULL, x3=NULL, x4=NULL, x5=NULL){
  #P(X1)
  if(!is.null(x1) & is.null(x2)){
    if(x1 == 1){
      return(3/4)
    } else {
      return(1/4)
    }
  }
  
  # P(X2 | X1)
  if (!is.null(x2) & !is.null(x1) & is.null(x3)){
    if (x2 == 1 & x1 == 1){
      return(4/5)
    } else if (x2 == 0 & x1 == 1){
      return(1/5)
    } else if (x2 == 1 & x1 == 0){
      return(2/3)
    } else if (x2 == 0 & x1 == 0){
      return(1/3)
    }
  }
  
  # P(X3 | X2)
  if (!is.null(x3) & !is.null(x2) & is.null(x4)){
    if (x3 == 1 & x2 == 1){
      return(5/7)
    } else if (x3 == 0 & x2 == 1){
      return(2/7)
    } else if (x3 == 1 & x2 == 0){
      return(1/3)
    } else if (x3 == 0 & x2 == 0){
      return(2/3)
    }
  }

  # P(X4 | X3)
  if (!is.null(x4) & !is.null(x3) & is.null(x5)){
    if (x4 == 1 & x3 == 1){
      return(3/5)
    } else if (x4 == 0 & x3 == 1){
      return(2/5)
    } else if (x4 == 1 & x3 == 0){
      return(2/5)
    } else if (x4 == 0 & x3 == 0){
      return(3/5)
    }
  }
  
  # P(X5 | X4)
  if (!is.null(x5) & !is.null(x4)){
    if (x5 == 1 & x4 == 1){
      return(1/2)
    } else if (x5 == 0 & x4 == 1){
      return(1/2)
    } else if (x5 == 1 & x4 == 0){
      return(7/9)
    } else if (x5 == 0 & x4 == 0){
      return(2/9)
    }
  }
}

# hardcoding the potential values
potentials = array(dim = c(2, 2, 4), dimnames 
                   = list(c("0", "1"), c("0", "1"), c("Psi12", "Psi23", "Psi34", "Psi45")))

potentials[1,1,1] = p(x1=0, x2=0) * p(x1=0)
potentials[1,2,1] = p(x1=0, x2=1) * p(x1=0)
potentials[2,1,1] = p(x1=1, x2=0) * p(x1=1)
potentials[2,2,1] = p(x1=1, x2=1) * p(x1=1)

potentials[1,1,2] = p(x2=0, x3=0)
potentials[1,2,2] = p(x2=0, x3=1)
potentials[2,1,2] = p(x2=1, x3=0)
potentials[2,2,2] = p(x2=1, x3=1)

potentials[1,1,3] = p(x3=0, x4=0)
potentials[1,2,3] = p(x3=0, x4=1)
potentials[2,1,3] = p(x3=1, x4=0)
potentials[2,2,3] = p(x3=1, x4=1)

potentials[1,1,4] = p(x4=0, x5=0)
potentials[1,2,4] = p(x4=0, x5=1)
potentials[2,1,4] = p(x4=1, x5=0)
potentials[2,2,4] = p(x4=1, x5=1)

potentials

# data structure to keep alpha messages
alpha_messages = array(dim = c(2,1,5), dimnames = list(c("0", "1"), c("mu_alpha"), c("x1", "x2", "x3", "x4", "x5")))

# initialiaing
alpha_messages[,,1] = c(1,1)

# computing the alpha message values being passed
for (i in 2:5){
  alpha_messages[,,i] = alpha_messages[,,i-1] %*% potentials[,,i-1]
}

alpha_messages

# data structure to keep beta messages
beta_messages = array(dim = c(2,1,5), dimnames = list(c("0", "1"), c("mu_beta"), c("x1", "x2", "x3", "x4", "x5")))

# initialiaing
beta_messages[,,5] = c(1,1)

# computing the alpha message values being passed
for (i in 4:1){
  beta_messages[,,i] = potentials[,,i] %*% beta_messages[,,i+1]
}

beta_messages

# multiplying forward and backward messages

marginal_probabilities = alpha_messages * beta_messages

z = sum(marginal_probabilities[,,3])
z
