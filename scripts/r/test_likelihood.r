# tree leaves (A, A) 
# newick: (A: 0.5, A: 0.5)

mu <- 1.0
t <- 0.5

delta <- (4.0 / 3.0) * t
p <- 0.25 * (1.0 + 3.0 * exp(-delta * mu))
q <- 0.25 * (1.0 - exp(-delta * mu))

P <- matrix(
  c(p, q, q, q,
    q, p, q, q,
    q, q, p, q,
    q, q, q, p), nrow=4, byrow=T)

err <- matrix(
  c(1.0 - 3.0 * ep, ep, ep, ep,
    ep, 1.0 - 3.0 * ep, ep, ep,
    ep, ep, 1.0 - 3.0 * ep, ep,
    ep, ep, ep, 1.0 - 3.0 * ep), nrow=4, byrow=T)

# P <- matrix(rep(q, 4 * 4), nrow=4)
# diag(P) <- rep(p, 4)
# err <- matrix(rep(ep, 4 * 4), nrow=4)
# diag(err) <- rep(1.0 - 3.0 * ep, 4)

p1 <- P %*% err[1,]
p2 <- P %*% err[2,]
p3 <- P %*% err[3,]
p4 <- P %*% err[4,]
prob <- 0.25 * (p1[1] ** 2 + p2[1] ** 2 + p3[1] ** 2 + p4[1] ** 2)
log(prob)
