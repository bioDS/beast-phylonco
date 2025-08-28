# Exponatiate matrix using several different methods:
#
# In Markov chain, the probability transition matrix P(t) is gained by:
# P(t) = e^(tQ)
# Numerically, there might be several problems to calculate e^(tQ), see:
# Nineteen Dubious Ways to Compute the Exponential of a Matrix, Twenty-Five Years Later (Moler & Van Loan, 2003)
#
# Note that this test suite is not perfect, it does not test every possible or extreme value and
# it does not make an attempt to understand the surface of solutions and possible problems.
# (e.g., the Hump in 19 Dubious Ways)
# All it does it that it tests matrix under several different parameter values with several
# different methods implemented in R. These results are then taken as a basic for BEAST test in Java. 
#
# Utilising packages "Matrix" and "expm"
# Making sure packages exist, but we will call particular functions directly
# due to namespace clashes
library("Matrix")
library("expm")
library("MASS")


set_up_matrix = function(k, a, b, g){
    data = c(
        0, 1, k, 1, 0, 0,
        1, 0, 1, k, a, 0,
        k, 1, 0, 1, 0, a,
        1, k, 1, 0, 0, 0,
        1, b, 1, k+g, 0, 0,
        k+g, 1, b, 1, 0, 0
        )
    Q = matrix(data, 6, 6, byrow=TRUE)
    for(i in 1:6){
        Q[i,i] = -sum(Q[i,])
        }
    Q   
    }


gsolve = function(A, b){
    ginv(A) %*% b
    }


stationary_distribution = function(Q){
    A = rbind(t(Q), 1)
    b = c(rep(0, nrow(Q)), 1)
    pi = c(gsolve(A, b))
    pi
    }


normalise = function(Q){
    pi = stationary_distribution(Q)
    mu = -sum(pi * diag(Q))
    Q/mu
    }


print.test_result = function(x, digits=10){
    x = x$results[[1]]
    x = formatC(x, digits=10, format="f")
    write.table(x, "", col.names=FALSE, row.names=FALSE, quote=FALSE, eol=",\n", sep=", ")
    }

expm = function(k, a, b, g, t, method){
    Q = set_up_matrix(k, a, b, g)
    Q = normalise(Q)
    expm_methods = c("Higham08.b", "AlMohy-Hi09", "Ward77", "PadeRBS") 
    methods = c("Matrix", expm_methods)

    if(!method %in% c("Matrix", expm_methods)){
        stop("Please use one of following methods: ", methods)
        }
    if(method == "Matrix"){
        Pt = as.matrix(Matrix::expm(t*Q))
        }
    if(method %in% expm_methods){
        Pt = expm::expm(t*Q, method=method, order=12)
        }
    return(Pt)
    }


test = function(k, a, b, g, t, tolerance=1e-10){
    results = list()
    methods = c("Matrix", "Higham08.b", "AlMohy-Hi09", "Ward77", "PadeRBS")
    for(method in methods){
        results[[method]] = expm(k, a, b, g, t, method)
        }

    # check if they are identical
    n = length(methods)
    comparison = matrix(NA, n, n)
    colnames(comparison) = methods
    rownames(comparison) = methods     
    for(i in 1:n){
        for(j in i:n){
            comparison[i, j] = all.equal(
                results[[i]], results[[j]],
                tolerance=tolerance, check.attributes=FALSE
                )
            comparison[j, i] = comparison[i, j]
            }
        }

    if(!all(comparison)){
        message("WARNING: possible numerical instability (see comparison)")
        }

    structure(
        list("results"=results, "comparison"=comparison),
        class = "test_result"
        )
    }


# Tests:
# k=a=b=g=t=1:
test1a = test(1,1,1,1,1)
# results:
#    0.4952638937, 0.1509017700, 0.1509017700, 0.1701551334, 0.0163887164, 0.0163887164,
#    0.1701551334, 0.3963929199, 0.1509017700, 0.1808145558, 0.0853469045, 0.0163887164,
#    0.1808145558, 0.1509017700, 0.3963929199, 0.1701551334, 0.0163887164, 0.0853469045,
#    0.1701551334, 0.1509017700, 0.1509017700, 0.4952638937, 0.0163887164, 0.0163887164,
#    0.1701551334, 0.1509017700, 0.1509017700, 0.2497727439, 0.2618798662, 0.0163887164,
#    0.2497727439, 0.1509017700, 0.1509017700, 0.1701551334, 0.0163887164, 0.2618798662,

#
# k=a=b=g=t, t=5
test1b = test(1,1,1,1,5)
# results:
#    0.2622804328, 0.1998216762, 0.1998216762, 0.2586484536, 0.0397138805, 0.0397138805,
#    0.2586484536, 0.2007132950, 0.1998216762, 0.2601365405, 0.0409661541, 0.0397138805,
#    0.2601365405, 0.1998216762, 0.2007132950, 0.2586484536, 0.0397138805, 0.0409661541,
#    0.2586484536, 0.1998216762, 0.1998216762, 0.2622804328, 0.0397138805, 0.0397138805,
#    0.2586484536, 0.1998216762, 0.1998216762, 0.2613888141, 0.0406054993, 0.0397138805,
#    0.2613888141, 0.1998216762, 0.1998216762, 0.2586484536, 0.0397138805, 0.0406054993,



# k=2, a=2, b=1, w=1, t=10
test2 = test(2,2,1,1,10)
# results:
#    0.2600943979, 0.1799686021, 0.1800318348, 0.2599056899, 0.0599684564, 0.0600310189,
#    0.2599056899, 0.1800308446, 0.1799686021, 0.2600940475, 0.0600323595, 0.0599684564,
#    0.2600940475, 0.1799686021, 0.1800308446, 0.2599056899, 0.0599684564, 0.0600323595,
#    0.2599056899, 0.1800318348, 0.1799686021, 0.2600943979, 0.0600310189, 0.0599684564,
#    0.2599056899, 0.1800311645, 0.1799686021, 0.2600947178, 0.0600313693, 0.0599684564,
#    0.2600947178, 0.1799686021, 0.1800311645, 0.2599056899, 0.0599684564, 0.0600313693,

# k=5, a=3, b=1, g=4, t=0.01
test3 = test(5,3,1,4, 0.01)
# results:
#    0.9919057846, 0.0011589132, 0.0057623071, 0.0011609356, 0.0000020162, 0.0000100434,
#    0.0011609356, 0.9884512168, 0.0011589132, 0.0057803839, 0.0034465344, 0.0000020162,
#    0.0057803839, 0.0011589132, 0.9884512168, 0.0011609356, 0.0000020162, 0.0034465344,
#    0.0011609356, 0.0057623071, 0.0011589132, 0.9919057846, 0.0000100434, 0.0000020162,
#    0.0011609356, 0.0011803191, 0.0011589132, 0.0103623719, 0.9861354441, 0.0000020162,
#    0.0103623719, 0.0011589132, 0.0011803191, 0.0011609356, 0.0000020162, 0.9861354441,



