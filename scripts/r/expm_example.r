# Example matrix exponentiation of ternary Q matrix from SiFit paper

library(expm)
op <- options(digits=12)
t <- 10
lambdaD <- 2
lambdaL <- 3
lambdaSum <- lambdaD + lambdaL
Q <- matrix(c(
	-1, 1, 0,
	lambdaSum / 2, -lambdaSum, lambdaSum / 2,
	0, lambdaD, -lambdaD
), nrow=3, byrow=TRUE)

expm(Q*t)


for (t in 1:10) {
  P <- expm(Q*t)
  print(P)
}