library('expm')

# translate from state index to genotype states
translate_seq <- function(seq, state) {
    states <- rep(state[1], length(seq))    
    for (i in 1:length(seq)) {
        states[i] = state[seq[i]]        
    }
    states
}

# generate random sequence
generate_seq <- function(l, state) {
    n_states <- length(state)
    seq <- runif(l, min=1, max=n_states+1)
    floor(seq)
}

# get the Q transition matrix
get_Q <- function(lambdaD, lambdaL) {
    lambdaSum <- lambdaD + lambdaL
    matrix(c(-1, 1, 0,
             lambdaSum/2, -lambdaSum, lambdaSum/2,
             0, lambdaD, -lambdaD), byrow=TRUE, ncol=3)
}

# get the P matrix P(t) = exp(Qt)
get_P <- function(Q, t) { 
    expm(Q * t)
}

# mutate single genotype
# g is the genotype index from the state set
# P is the matrix P(t) = exp(Qt)
mutate <- function(g, P) {
    prob <- P[g, ]
    cumprob <- cumsum(prob)
    r <- runif(1)
    for (i in 1:length(prob)) {
        if (r < cumprob[i]) {
            return (i)
        }
    }
}

# mutate sequence
mutate_seq <- function(seq, P) {
    child_seq <- seq
    for (i in 1:length(seq)) {
        child_seq[i] <- mutate(child_seq[i], P)
    }
    child_seq
}

mutate_tree <- function(tree, Q) {
    root <- get_root()
    children <- root.get_children()        
    branch <- parent.get_time() - child.get_time()
    P <- get_P(Q, branch)
    mutate_seq(parent.get_seq(), P)
}

# example 1: generate 1 child
example_1 <- function() {
    lambdaD <- 3
    lambdaL <- 2
    l <- 20
    t <- 5

    state = c(0, 1, 2)
    Q <- get_Q(lambdaD, lambdaL)
    P <- get_P(Q, t)

    ancestor <- rep(1, l)
    child <- mutate_seq(ancestor, P)

    cat("example 1 - generate 1 child\n")
    cat("ancestor sequence:\n")
    cat(translate_seq(ancestor, state))
    cat("\nchild sequence:\n")
    cat(translate_seq(child, state))
}

# example 2: generate 50 children
example_2 <- function() {
    lambdaD <- 3
    lambdaL <- 2
    l <- 20
    t <- 5

    state = c(0, 1, 2)
    Q <- get_Q(lambdaD, lambdaL)
    P <- get_P(Q, t)

    ancestor <- rep(1, l) 

    cat("\n\nexample 2 - generate 50 children\n")
    cat("ancestor sequence:\n")
    cat(translate_seq(ancestor, state))
    cat("\nchild sequence:\n")

    for (i in 1:50) {
        child <- mutate_seq(ancestor, P)
        cat(translate_seq(child, state))
        cat("\n")
    }
}