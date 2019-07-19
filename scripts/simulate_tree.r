source("simulate_seq.r")

# simulate tree with newick format
#  (((A:t2, B:t2)D:t1, C:t1 + t2)E);
# with sequence length l = 20
# and substitution model parameters
#  lambdaD = 3 and lambdaL = 2
# the default tree is 
#  (((A:0.5, B:0.5)D:0.5, C:1)E);
simulate_tree_small <- function(lambdaD=3, lambdaL=2, l=20, t1=0.5, t2=0.5) {    
    state = c("0", "1", "2")
    Q <- get_Q(lambdaD, lambdaL)
    
    e <- rep(1, l)
    d <- mutate_seq(e, get_P(Q, t1))
    a <- mutate_seq(d, get_P(Q, t2))
    b <- mutate_seq(d, get_P(Q, t2))
    c <- mutate_seq(e, get_P(Q, t1 + t2))
    
    a_seq <- translate_seq(a, state)
    b_seq <- translate_seq(b, state)
    c_seq <- translate_seq(c, state)
    d_seq <- translate_seq(d, state)
    e_seq <- translate_seq(e, state)
    
    list(a_seq, b_seq, c_seq)
}

# generate 100 trees with fixed lambdas
#  lambdaD = 3
#  lambdaL = 2
sim_1 <- function(filename, N=100, lambdaD=3, lambdaL=2, seed=777) {
    set.seed(seed)
    cat("writing sim_1() to file: ", filename, "\n")
    sink(filename)
    cat("tree, lambdaD, lambdaL, node, sequence\n")
    for (i in 1:N) {
        tree <- simulate_tree_small()
        cat(i, lambdaD, lambdaL, "a", "", sep=", ")
        cat(tree[[1]], "\n", sep="")
        cat(i, lambdaD, lambdaL, "b", "", sep=", ")
        cat(tree[[2]], "\n", sep="")
        cat(i, lambdaD, lambdaL, "c", "", sep=", ")
        cat(tree[[3]], "\n", sep="")
    }
    sink()
}

# generate 100 trees with lambdas
#  lambdaD ~ LogNorm(mu = -1)
#  lambdaL ~ LogNorm(mu = 0.5)
sim_2 <- function(filename, N=100, seed=777) {
    set.seed(seed)
    cat("writing sim_2() to file: ", filename, "\n")
    sink(filename)
    cat("tree, lambdaD, lambdaL, node, sequence\n")
    for (i in 1:N) {
        lambdaD <- rlnorm(1, -1)
        lambdaL <- rlnorm(1, 0.5)
        tree <- simulate_tree_small(lambdaD, lambdaL)
        cat(i, lambdaD, lambdaL, "a", "", sep=", ")
        cat(tree[[1]], "\n", sep="")
        cat(i, lambdaD, lambdaL, "b", "", sep=", ")
        cat(tree[[2]], "\n", sep="")
        cat(i, lambdaD, lambdaL, "c", "", sep=", ")
        cat(tree[[3]], "\n", sep="")
    }
    sink()
}

# run simulations
sim_1("tree_small_fixed.csv")
sim_2("tree_small_lognorm.csv")