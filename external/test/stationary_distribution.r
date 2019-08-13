# Calculate stationary distribution:
# pi Q = 0
# where pi is stationary distribution and Q is rate matrix.
#
# In BEAST, stationary distribution is used either:
# 1) To extend rate matrix:
#   -- Q = R*F where F are frequencies of stationary distribution
#   -- This makes stationary distribution property of data instead of rate matrix R
#   -- not used in this case due to non time-reversibility
# 2) To normalize rate matrix
#   -- rate matrix Q is normalized to a single expected transition per time
#   -- Q = mu Q'
# 3) As nucleotide frequencies near the root
#
# However, I cannot get the solution symbolically, both SymPy and Yacas fails to deliver an answer.
# Numerical solution however works fine
#
# The numerical solution is calculated using an generalised inverse of transposed A.
#
# The methylationHKY matrix is:
#
#    A C G T C' G'
# A  - 1 k 1 0  0
# C  1 - 1 k a  0
# T  k 1 - 1 0  a
# G  1 k 1 - 0  0
# C' 1 b 1 k+g - 0
# G' k+g 1 b 1 0 -
#
# where - is negative sum of other elements in a row
#
library(MASS)


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

stationary = function(k, a, b, g){
    Q = set_up_matrix(k, a, b, g)
    pi = stationary_distribution(Q)
    pi
    }

# Some examples:
# k=a=b=g=1, tested by hand
# result: 0.26 0.20 0.20 0.26 0.04 0.04
test1 = stationary(1, 1, 1, 1)

# k=2, a=2, b=1, w=1, also tested by hand
# result: 0.26 0.18 0.18 0.26 0.06 0.06
test2 = stationary(2, 2, 1, 1)

# k=5, a=3, b=1, g=4
# result: 0.265625 0.187500 0.187500 0.265625 0.046875 0.046875
test3 = stationary(5, 3, 1, 4)
