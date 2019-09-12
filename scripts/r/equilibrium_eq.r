library(matlib)

# Modified SiFit3 matrix
rows <- 3
lambdaD <- 3
lambdaL <- 2
lambdaSum <- lambdaD + lambdaL
Q <- matrix(c(
     -1, 1, 0,
     1/6 + lambdaSum / 2, -lambdaSum - 1/6 - 1/2, 1/2 + lambdaSum / 2,
     0, lambdaD, -lambdaD
 ), nrow=rows, byrow=T)

A <- rbind(t(Q) , rep(1, rows))
b <- c(rep(0, rows), 1)
Solve(A, b)

# K80 diploid matrix
kappa <- 3
lambdaL <- 2
rows <- 10

Q <- matrix(c(
	0, 1, kappa, 1, 0, 0, 0, 0, 0, 0,
	1, 0, 1 / 2, kappa / 2, 1, kappa / 2, 1 / 2, 0, 0, 0,
	kappa, 1 / 2, 0, 1 / 2, 0, 1 / 2, 0, kappa, 1 / 2, 0,
	1, kappa / 2, 1 / 2, 0, 0, 0, 1 / 2, 0, kappa / 2, 1,
	0, 1, 0, 0, 0, 1, kappa, 0, 0, 0,
	0, kappa / 2, 1 / 2, 0, 1, 0, 1 / 2, 1, kappa / 2, 0,
	0, 1 / 2, 0, 1 / 2, kappa, 1 / 2, 0, 0, 1 / 2, kappa,
	0, 0, kappa, 0, 0, 1, 0, 0, 1, 0,
	0, 0, 1 / 2, kappa / 2, 0, kappa / 2, 1 / 2, 1, 0, 1,
	0, 0, 0, 1, 0, 0, kappa, 0, 1, 0
 ), nrow=rows, byrow=T)

# loss of one allele
Q[2,1] <- Q[2,1] + lambdaL
Q[3,1] <- Q[3,1] + lambdaL
Q[4,1] <- Q[4,1] + lambdaL
Q[2,5] <- Q[2,5] + lambdaL
Q[6,5] <- Q[6,5] + lambdaL
Q[7,5] <- Q[7,5] + lambdaL
Q[3,8] <- Q[3,8] + lambdaL
Q[6,8] <- Q[6,8] + lambdaL
Q[9,8] <- Q[9,8] + lambdaL
Q[4,10] <- Q[4,10] + lambdaL
Q[7,10] <- Q[7,10] + lambdaL
Q[9,10] <- Q[9,10] + lambdaL

# diagonals
for (i in 1:rows) {
	Q[i,i] <- -sum(Q[i,]) + Q[i,i]
}

A <- rbind(t(Q) , rep(1, rows))
b <- c(rep(0, rows), 1)
Solve(A, b)
