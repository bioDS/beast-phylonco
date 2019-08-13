import numpy
import scipy.linalg





def set_up_matrix(k, a, b, g):
    mat = numpy.matrix([
        [0, 1, k, 1, 0, 0],
        [1, 0, 1, k, a, 0],
        [k, 1, 0, 1, 0, a],
        [1, k, 1, 0, 0, 0],
        [1, b, 1, k+g, 0, 0],
        [k+g, 1, b, 1, 0, 0]
        ])
    for i in range(mat.shape[0]):
        mat[i, i] = -mat[i].sum()
    return mat


def stationary_distribution(Q):
    A = numpy.concatenate(
        (Q.T, numpy.matrix(numpy.ones(Q.shape[0])))
        )
    b = numpy.matrix(
        numpy.append(numpy.zeros(Q.shape[0]), 1)
        )
    pi = gsolve(A,b.T)
    return pi


def gsolve(A, b):
    return numpy.matmul(numpy.linalg.pinv(A), b) 


def normalise(Q):
    pi = stationary_distribution(Q)
    mu = -numpy.multiply(pi, Q.diagonal().T).sum()
    return Q / mu


def expm(k, a, b, g, t):
    Q = set_up_matrix(k, a, b, g)
    Q = normalise(Q)
    Pt = scipy.linalg.expm(Q*t)
    return Pt


def main():
    Pt = expm(1,1,1,1,1)
    print(Pt)


if __name__ == "__main__":
    main()
