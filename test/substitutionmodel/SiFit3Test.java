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
     *	 -1, 1, 0,
     *	 lambdaSum / 2, -lambdaSum, lambdaSum / 2,
     *	 0, lambdaD, -lambdaD
     * ), nrow=3, byrow=TRUE)
     * expm(Q*t)
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
                0.526315996407, 0.210526278501, 0.263157725092,
                0.526315696252, 0.210526332588, 0.263157971160,
                0.526315450185, 0.210526376928, 0.263158172887
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
     *	 -1, 1, 0,
     *	 lambdaSum / 2, -lambdaSum, lambdaSum / 2,
     *	 0, lambdaD, -lambdaD
     * ), nrow=3, byrow=TRUE)
     * expm(Q*t)
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
                    0.9148390890826, 0.0754940428491, 0.00966686806831,
                    0.1887351071226, 0.6321966538230, 0.17906823905434,
                    0.0193337361366, 0.1432545912435, 0.83741167261991
            };
        }
    };

    Instance[] allTests = {test0, test1};

    public void  testSiFit3() throws Exception {
        for (Instance test: allTests) {
            SiFit3 siFit = new SiFit3();
            siFit.initByName(
                    "lambdaD", new RealParameter(test.getLambdaD().toString()),
                    "lambdaL", new RealParameter(test.getLambdaL().toString())
            );

            double distance = test.getDistance();
            double[] expectedResult = test.getExpectedResult();
            double[] matrix = new double[3 * 3];

            siFit.getTransitionProbabilities(null, distance, 0, 1, matrix);
            for (int k = 0; k < matrix.length; ++k) {
                assertEquals(matrix[k], expectedResult[k], 1e-10);
                System.out.println(k + " : " + (matrix[k] - expectedResult[k]));
                System.out.println(expectedResult[k]);
                System.out.println(matrix[k] + "\n");
            }
        }
    }

}
