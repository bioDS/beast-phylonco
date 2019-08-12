package substitutionmodel;

import beast.core.Description;
import beast.core.parameter.RealParameter;
import junit.framework.TestCase;

@Description("Test SiFit3 matrix exponentiation")
public class SiFit3Test extends TestCase {

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
     *     -1, 1, 0,
     *     lambdaSum / 2, -lambdaSum, lambdaSum / 2,
     *     0, lambdaD, -lambdaD
     * ), nrow=3, byrow=T)
     * x <- 1 + 2 / lambdaSum + 1 / lambdaD
     * pi0 <- 1 / x
     * pi1 <- 2 / (x * lambdaSum)
     * pi2 <- 1 / (x * lambdaD)
     * freq <- c(pi0, pi1, pi2)
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
                0.526735561965, 0.210450674425, 0.262813763610,
                0.526126686062, 0.210560391486, 0.263312922452,
                0.525627527220, 0.210650337962, 0.263722134819
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
     *     -1, 1, 0,
     *     lambdaSum / 2, -lambdaSum, lambdaSum / 2,
     *     0, lambdaD, -lambdaD
     * ), nrow=3, byrow=T)
     * x <- 1 + 2 / lambdaSum + 1 / lambdaD
     * pi0 <- 1 / x
     * pi1 <- 2 / (x * lambdaSum)
     * pi2 <- 1 / (x * lambdaD)
     * freq <- c(pi0, pi1, pi2)
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
                    0.95614090897870, 0.0413688413513, 0.00249024967001,
                    0.10342210337822, 0.7956460429136, 0.10093185370821,
                    0.00498049934002, 0.0807454829666, 0.91427401769341
            };
        }
    };

    Instance[] allTests = {test0, test1};

    public void  testSiFit3() throws Exception {
        for (Instance test: allTests) {
            SiFit3 model = new SiFit3();
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
