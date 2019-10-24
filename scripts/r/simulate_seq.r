# requirements
require('expm')

# Implementation of the SiFit substitution model for ternary data 
# References: Zafar et al (2017).
# ---------------------------------------------------------------

states <- c("0", "1", "2")

# get the equilibrium state frequencies
getPi <- function(lambdaD, lambdaL) {
  lambdaSum <- lambdaD + lambdaL
  x <- 1 + 2 / lambdaSum + 1 / lambdaD
  pi0 <- 1 / x
  pi1 <- 2 / (x * lambdaSum) 
  pi2 <- 1 / (x * lambdaD)
  c(pi0, pi1, pi2)
}

# get the normalization factor for the Q matrix
getBeta <- function(lambdaD, lambdaL) {
  freq <- getPi(lambdaD, lambdaL)
  diag <- c(1, lambdaD + lambdaL, lambdaD)
  as.vector(1 / (freq %*% diag))
}

# get the normalized Q transition matrix
getQ <- function(lambdaD, lambdaL) {
  lambdaSum <- lambdaD + lambdaL
  Q <- matrix(c(-1, 1, 0,
    lambdaSum/2, -lambdaSum, lambdaSum/2,
    0, lambdaD, -lambdaD), byrow=T, ncol=3)
  beta <- getBeta(lambdaD, lambdaL)
  beta * Q
}

# get the P matrix P(t) = exp(Qt)
getP <- function(Q, t) { 
  expm(Q * t)
}

# Simulation of sequences
# ---------------------------------------------------------------

# sample from a given probability vector
sampleFrom <- function(prob) {  
  cumProb <- cumsum(prob)
  r <- runif(1)
  for (i in 1:length(cumProb)) {
    if (r < cumProb[i]) {
      return (i)
    }
  }
}

# generate sequence from frequencies
generateSeq <- function(l, freq) {
  seq <- rep(0, l)
  for (i in 1:l) {
    seq[i] = sampleFrom(freq)
  }
  seq
}

# mutate sequence
mutateSeq <- function(seq, P) {
  childSeq <- seq
  for (i in 1:length(seq)) {
    prob <- P[childSeq[i], ]
    childSeq[i] <- sampleFrom(prob)
  }
  childSeq
}

# translate from state index to actual states
translateSeq <- function(seq) {
  stateSeq <- rep(states[1], length(seq))  
  for (i in 1:length(seq)) {
    stateSeq[i] = states[seq[i]]    
  }
  stateSeq
}