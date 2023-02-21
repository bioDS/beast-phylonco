package phylonco.beast.evolution.substitutionmodel;

import beast.base.core.Description;
import beast.base.inference.parameter.RealParameter;
import org.junit.Test;
import phylonco.beast.evolution.substitutionmodel.BinarySubstitutionModel;

import static org.junit.Assert.assertArrayEquals;

@Description("Test Binary substitution model matrix exponentiation")
public class BinarySubstitutionModelTest {

    private static double DELTA = 1e-10;

    private BinarySubstitutionModel model;
    private int nrOfStates;

    public void setupModel(Double lambda) {
        model = new BinarySubstitutionModel();
        model.initByName(
                "lambda", new RealParameter(lambda.toString())
        );
        nrOfStates = model.getStateCount();
    }

    /**
     * results obtained from running the following code in R:
     *
     * library(expm)
     * op <- options(digits=12)
     * t <- 10
     * lambda <- 0.3
     * Q <- matrix(c(
     *     -1, 1,
     *     lambda, -lambda
     * ), nrow=2, byrow=T)
     * pi0 <- lambda / (lambda + 1)
     * pi1 <- 1 / (lambda + 1)
     * freq <- c(pi0, pi1)
     * diag <- -diag(Q)
     * beta <- as.vector(1 / (freq %*% diag))
     * expm(beta * Q * t)
     *
     */
    @Test
    public void testLambdaSmall() {
        double t = 10;
        setupModel(0.3);
        double[] expected = new double[] {
                0.230769230770, 0.769230769230,
                0.230769230769, 0.769230769231
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
     * lambda <- 3
     * Q <- matrix(c(
     *     -1, 1,
     *     lambda, -lambda
     * ), nrow=2, byrow=T)
     * pi0 <- lambda / (lambda + 1)
     * pi1 <- 1 / (lambda + 1)
     * freq <- c(pi0, pi1)
     * diag <- -diag(Q)
     * beta <- as.vector(1 / (freq %*% diag))
     * expm(beta * Q * t)
     *
     */
    @Test
    public void testlambdaLarge() {
        double t = 10;
        setupModel(3.0);
        double[] expected = new double[] {
                0.750000000001, 0.249999999999,
                0.749999999998, 0.250000000002
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
     * lambda <- 2
     * Q <- matrix(c(
     *     -1, 1,
     *     lambda, -lambda
     * ), nrow=2, byrow=T)
     * pi0 <- lambda / (lambda + 1)
     * pi1 <- 1 / (lambda + 1)
     * freq <- c(pi0, pi1)
     * diag <- -diag(Q)
     * beta <- as.vector(1 / (freq %*% diag))
     * expm(beta * Q * t)
     *
     */
    @Test
    public void testTransitionLong() {
        double t = 10;
        setupModel(2.0);
        double[] expected = new double[] {
                0.666666666723, 0.333333333277,
                0.666666666554, 0.333333333446
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
     * lambda <- 2
     * Q <- matrix(c(
     *     -1, 1,
     *     lambda, -lambda
     * ), nrow=2, byrow=T)
     * pi0 <- lambda / (lambda + 1)
     * pi1 <- 1 / (lambda + 1)
     * freq <- c(pi0, pi1)
     * diag <- -diag(Q)
     * beta <- as.vector(1 / (freq %*% diag))
     * expm(beta * Q * t)
     *
     */
    @Test
    public void testTransitionShort() {
        double t = 0.1;
        setupModel(2.0);
        double[] expected = new double[] {
                0.932838739586, 0.0671612604135,
                0.134322520827, 0.8656774791729
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
     * lambda <- 2
     * Q <- matrix(c(
     *     -1, 1,
     *     lambda, -lambda
     * ), nrow=2, byrow=T)
     * pi0 <- lambda / (lambda + 1)
     * pi1 <- 1 / (lambda + 1)
     * freq <- c(pi0, pi1)
     * diag <- -diag(Q)
     * beta <- as.vector(1 / (freq %*% diag))
     * result <- expm(beta * Q * t)
     * result[1, ]
     *
     */
    @Test
    public void testEquilibrium() {
        setupModel(2.0);
        double[] expected = new double[] {
                0.666666666723, 0.333333333277
        };
        double[] observed = model.getFrequencies();
        assertArrayEquals(expected, observed, DELTA);
    }
}
