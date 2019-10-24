package substitutionmodel;

import beast.core.Description;
import beast.core.parameter.RealParameter;

import org.junit.*;
import static org.junit.Assert.*;

@Description("Test SiFit2 matrix exponentiation")
public class SiFit2Test {

    private static double DELTA = 1e-10;

    private SiFit2 model;
    private int nrOfStates;

    public void setupModel(Double lambdaD, Double lambdaL) {
        model = new SiFit2();
        model.initByName(
                "lambdaD", new RealParameter(lambdaD.toString()),
                "lambdaL", new RealParameter(lambdaL.toString())
        );
        nrOfStates = model.getStateCount();
    }

    /**
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=12)
     * t <- 10
     * lambdaD <- 2
     * lambdaL <- 3
     * lambdaSum <- lambdaD + lambdaL
     * Q <- matrix(c(
     *     -1, 1,
     *     lambdaSum / 2, -lambdaSum / 2
     * ), nrow=2, byrow=T)
     * x <- 1 + lambdaSum / 2
     * pi0 <- lambdaSum / (2 * x)
     * pi1 <- 1 / x
     * freq <- c(pi0, pi1)
     * diag <- -diag(Q)
     * beta <- as.vector(1 / (freq %*% diag))
     * expm(beta * Q * t)
     *
     */
    @Test
    public void testTransitionLong() {
        double t = 10;
        setupModel(2.0, 3.0);
        double[] expected = new double[] {
                0.714285714292, 0.285714285708,
                0.714285714269, 0.285714285731
        };
        double[] observed = new double[nrOfStates * nrOfStates];
        model.getTransitionProbabilities(null, t, 0, 1, observed);
        assertArrayEquals(expected, observed, DELTA);
    }

    /**
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=12)
     * t <- 0.1
     * lambdaD <- 2
     * lambdaL <- 3
     * lambdaSum <- lambdaD + lambdaL
     * Q <- matrix(c(
     *     -1, 1,
     *     lambdaSum / 2, -lambdaSum / 2
     * ), nrow=2, byrow=T)
     * x <- 1 + lambdaSum / 2
     * pi0 <- lambdaSum / (2 * x)
     * pi1 <- 1 / x
     * freq <- c(pi0, pi1)
     * diag <- -diag(Q)
     * beta <- as.vector(1 / (freq %*% diag))
     * expm(beta * Q * t)
     *
     */
    @Test
    public void testTransitionShort() {
        double t = 0.1;
        setupModel(2.0, 3.0);
        double[] expected = new double[] {
                0.937915582355, 0.0620844176452,
                0.155211044113, 0.8447889558870
        };
        double[] observed = new double[nrOfStates * nrOfStates];
        model.getTransitionProbabilities(null, t, 0, 1, observed);
        assertArrayEquals(expected, observed, DELTA);
    }

    /**
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=12)
     * t <- 10
     * lambdaD <- 2
     * lambdaL <- 3
     * lambdaSum <- lambdaD + lambdaL
     * Q <- matrix(c(
     *     -1, 1,
     *     lambdaSum / 2, -lambdaSum / 2
     * ), nrow=2, byrow=T)
     * x <- 1 + lambdaSum / 2
     * pi0 <- lambdaSum / (2 * x)
     * pi1 <- 1 / x
     * freq <- c(pi0, pi1)
     * diag <- -diag(Q)
     * beta <- as.vector(1 / (freq %*% diag))
     * result <- expm(beta * Q * t)
     * result[1, ]
     *
     */
    @Test
    public void testEquilibrium() {
        setupModel(2.0, 3.0);
        double[] expected = new double[] {
                0.714285714292, 0.285714285708
        };
        double[] observed = model.getFrequencies();
        assertArrayEquals(expected, observed, DELTA);
    }
}