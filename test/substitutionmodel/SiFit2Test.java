package substitutionmodel;

import beast.core.Description;
import beast.core.parameter.RealParameter;
import junit.framework.TestCase;

@Description("Test SiFit2 matrix exponentiation")
public class SiFit2Test extends TestCase {

    public interface Instance {
        Double getLambdaD();
        Double getLambdaL();

        double getDistance();
        double[] getExpectedResult();
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
    protected Instance test0 = new Instance() {
        @Override
        public Double getLambdaD() {
            return 2.0;
        }

        @Override
        public Double getLambdaL() {
            return 3.0;
        }

        @Override
        public double getDistance() {
            return 10;
        }

        @Override
        public double[] getExpectedResult() {
            return new double[] {
                    0.714285714292, 0.285714285708,
                    0.714285714269, 0.285714285731
            };
        }
    };

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
    protected Instance test1 = new Instance() {
        @Override
        public Double getLambdaD() {
            return 2.0;
        }

        @Override
        public Double getLambdaL() {
            return 3.0;
        }

        @Override
        public double getDistance() {
            return 0.1;
        }

        @Override
        public double[] getExpectedResult() {
            return new double[] {
                    0.937915582355, 0.0620844176452,
                    0.155211044113, 0.8447889558870
            };
        }
    };

    /**
     * expect probability to reach equilibrium frequencies when time t is large
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
    protected Instance test2 = new Instance() {
        @Override
        public Double getLambdaD() {
            return 2.0;
        }

        @Override
        public Double getLambdaL() {
            return 3.0;
        }

        @Override
        public double getDistance() {
            return 10;
        }

        @Override
        public double[] getExpectedResult() {
            SiFit2 model = new SiFit2();
            model.initByName(
                    "lambdaD", new RealParameter(getLambdaD().toString()),
                    "lambdaL", new RealParameter(getLambdaL().toString())
            );
            // expect equilibrium frequencies
            int nrStates = model.getStateCount();
            double[] freq = model.getFrequencies();
            System.out.println(freq);
            double[] result = new double[nrStates * nrStates];
            for (int i = 0; i < result.length; i++) {
                result[i] = freq[i % nrStates];
            }
            return result;
        }
    };

    Instance[] allTests = {test0, test1, test2};

    public void testSiFit2() throws Exception {
        for (Instance test: allTests) {
            SiFit2 model = new SiFit2();
            model.initByName(
                    "lambdaD", new RealParameter(test.getLambdaD().toString()),
                    "lambdaL", new RealParameter(test.getLambdaL().toString())
            );

            int nrStates = model.getStateCount();
            double distance = test.getDistance();
            double[] expectedResult = test.getExpectedResult();
            double[] matrix = new double[nrStates * nrStates];

            model.getTransitionProbabilities(null, distance, 0, 1, matrix);
            for (int k = 0; k < matrix.length; ++k) {
                assertEquals(matrix[k], expectedResult[k], 1e-10);
                System.out.println(k + " : " + (matrix[k] - expectedResult[k]));
            }
        }
    }
}
