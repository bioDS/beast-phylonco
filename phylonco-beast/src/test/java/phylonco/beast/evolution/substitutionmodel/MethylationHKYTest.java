package phylonco.beast.evolution.substitutionmodel;

import beast.base.core.Description;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import junit.framework.TestCase;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import phylonco.beast.TestUtils;

import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.MalformedURLException;

@Description("Test stationary distribution and matrix exponentiation for MethylationHKY matrix")
public class MethylationHKYTest extends TestCase {

    /** Helper function for printing matrices in array form
     *
     * The matrix stored in array is assumed to be in following form:
     * 0 1
     * 2 3
     * where the number is index of array, i.e.,:
     * array[0, 1, 2, 3]
     *
     *
     * @param matrix -- the matrix (in the form of unstructured arrays)
     * @param rows -- rows of matrix
     * @param columns -- columns of matrix
     */
    private void printMatrix(double[] matrix, int rows, int columns){
        for(int i=0; i < rows; i++){
            for(int j=0; j < columns; j++){
                System.out.print(matrix[i * rows + j] + ", ");
            }
            System.out.println();
        }
    }


    private void printMatrix(double[][] matrix, int rows, int columns){
        for(int i=0; i < matrix.length; i++){
            for(int j=0; j< columns; j++){
                System.out.print(matrix[i][j] + ", ");
            }
            System.out.println();
        }
    }

    class testData {
        RealParameter kappa;
        RealParameter alpha;
        RealParameter beta;
        RealParameter gamma;
        double[] result;
        double time;

        public testData(double[] input, double time, double[] result) {
            this.kappa = new RealParameter(input[0] + "");
            this.alpha = new RealParameter(input[1] + "");
            this.beta = new RealParameter(input[2] + "");
            this.gamma = new RealParameter(input[3] + "");
            this.time = time;
            this.result = result;
        }

        public testData(double[] input, double[] result){
            this(input, 0, result);
        }
    }

    @BeforeClass
    public static void setUpClass() {
        TestUtils.loadServices();
    }

    @Test
    public void testStationaryDistribution() {
        /* calculated numerically with R using Mooee-Penrose pseudoinversion, see:
         * external/test/stationary_distribution.r
         * for code and results
         */
        testData[] data = new testData[]{
                new testData(
                        new double[]{1, 1, 1, 1}, new double[]{0.26, 0.20, 0.20, 0.26, 0.04, 0.04}
                        ),
                new testData(
                        new double[]{2, 2, 1, 1}, new double[]{0.26, 0.18, 0.18, 0.26, 0.06, 0.06}
                        ),
                new testData(
                        new double[]{5, 3, 1, 4}, new double[]{0.265625, 0.187500, 0.187500, 0.265625, 0.046875, 0.046875}
                        )
                };

        for(testData test : data){
            MethylationHKY substitutionModel = new MethylationHKY();
            substitutionModel.initByName(
                    "kappa", test.kappa,
                    "alpha", test.alpha,
                    "beta", test.beta,
                    "gamma", test.gamma
                    );
            Assert.assertArrayEquals(test.result, substitutionModel.getFrequencies(), 1e-10);
        }
    }


    @Test
    public void testMatrixExponentiation() {
        /* calculated using R's Matrix and expm packages for several methods, see:
         * external/test/matrix_exponentiation.r
         */
        testData[] data = new testData[]{
            new testData(
                new double[]{1, 1, 1, 1}, 1, new double[]{
                    0.4952638937, 0.1509017700, 0.1509017700, 0.1701551334, 0.0163887164, 0.0163887164,
                    0.1701551334, 0.3963929199, 0.1509017700, 0.1808145558, 0.0853469045, 0.0163887164,
                    0.1808145558, 0.1509017700, 0.3963929199, 0.1701551334, 0.0163887164, 0.0853469045,
                    0.1701551334, 0.1509017700, 0.1509017700, 0.4952638937, 0.0163887164, 0.0163887164,
                    0.1701551334, 0.1509017700, 0.1509017700, 0.2497727439, 0.2618798662, 0.0163887164,
                    0.2497727439, 0.1509017700, 0.1509017700, 0.1701551334, 0.0163887164, 0.2618798662,

                }
            ),
            new testData(
                new double[]{1, 1, 1, 1}, 5, new double[]{
                    0.2622804328, 0.1998216762, 0.1998216762, 0.2586484536, 0.0397138805, 0.0397138805,
                    0.2586484536, 0.2007132950, 0.1998216762, 0.2601365405, 0.0409661541, 0.0397138805,
                    0.2601365405, 0.1998216762, 0.2007132950, 0.2586484536, 0.0397138805, 0.0409661541,
                    0.2586484536, 0.1998216762, 0.1998216762, 0.2622804328, 0.0397138805, 0.0397138805,
                    0.2586484536, 0.1998216762, 0.1998216762, 0.2613888141, 0.0406054993, 0.0397138805,
                    0.2613888141, 0.1998216762, 0.1998216762, 0.2586484536, 0.0397138805, 0.0406054993,
                }
            ),
            new testData(
                new double[]{2, 2, 1, 1}, 10, new double[]{
                    0.2600943979, 0.1799686021, 0.1800318348, 0.2599056899, 0.0599684564, 0.0600310189,
                    0.2599056899, 0.1800308446, 0.1799686021, 0.2600940475, 0.0600323595, 0.0599684564,
                    0.2600940475, 0.1799686021, 0.1800308446, 0.2599056899, 0.0599684564, 0.0600323595,
                    0.2599056899, 0.1800318348, 0.1799686021, 0.2600943979, 0.0600310189, 0.0599684564,
                    0.2599056899, 0.1800311645, 0.1799686021, 0.2600947178, 0.0600313693, 0.0599684564,
                    0.2600947178, 0.1799686021, 0.1800311645, 0.2599056899, 0.0599684564, 0.0600313693,

                }
            ),
            new testData(
                new double[]{5, 3, 1, 4}, 0.01, new double[]{
                    0.9919057846, 0.0011589132, 0.0057623071, 0.0011609356, 0.0000020162, 0.0000100434,
                    0.0011609356, 0.9884512168, 0.0011589132, 0.0057803839, 0.0034465344, 0.0000020162,
                    0.0057803839, 0.0011589132, 0.9884512168, 0.0011609356, 0.0000020162, 0.0034465344,
                    0.0011609356, 0.0057623071, 0.0011589132, 0.9919057846, 0.0000100434, 0.0000020162,
                    0.0011609356, 0.0011803191, 0.0011589132, 0.0103623719, 0.9861354441, 0.0000020162,
                    0.0103623719, 0.0011589132, 0.0011803191, 0.0011609356, 0.0000020162, 0.9861354441,
                }
            )
        };
        for (testData test : data){
            MethylationHKY substitutionModel = new MethylationHKY();
            substitutionModel.initByName(
                "kappa", test.kappa,
                    "alpha", test.alpha,
                    "beta", test.beta,
                    "gamma", test.gamma
            );
            double[] matrix = new double[6 * 6];
            substitutionModel.getTransitionProbabilities(null, test.time, 0, 1, matrix);
            //printMatrix(matrix, 6, 6);
            Assert.assertArrayEquals(test.result, matrix, 1e-9);
        }
    }



    @Test
    @SuppressWarnings("unchecked")
    public void testUnnormalizedRateMatrix() throws NoSuchMethodException, InvocationTargetException, IllegalAccessException {
        testData[] data = new testData[]{
            new testData(
                new double[]{1,1,1,1}, new double[]{
                    -3, 1, 1, 1, 0, 0,
                    1, -4, 1, 1, 1, 0,
                    1, 1, -4, 1, 0, 1,
                    1, 1, 1, -3, 0, 0,
                    1, 1, 1, 2, -5, 0,
                    2, 1, 1, 1, 0, -5,
                }
            ),
            new testData(
                new double[]{4, 3, 2, 1}, new double[]{
                    -6, 1, 4, 1, 0, 0,
                    1, -9, 1, 4, 3, 0,
                    4, 1, -9, 1, 0, 3,
                    1, 4, 1, -6, 0, 0,
                    1, 2, 1, 5, -9, 0,
                    5, 1, 2, 1, 0, -9,
                }
            )
        };
        for(testData test : data) {
            MethylationHKY substitutionModel = new MethylationHKY();
            substitutionModel.initByName(
                    "kappa", test.kappa,
                    "alpha", test.alpha,
                    "beta", test.beta,
                    "gamma", test.gamma
            );
            Method method = MethylationHKY.class.getDeclaredMethod("setupUnnormalizedRateMatrix");
            method.setAccessible(true);
            double[][] unnormalizedMatrix = (double[][]) method.invoke(substitutionModel);
            double[] flatMatrix = flattenMatrix(unnormalizedMatrix, 6, 6);
            //printMatrix(flatMatrix, 6, 6);
            Assert.assertArrayEquals(test.result, flatMatrix, 1e-10);
        }
    }


    @Test
    @SuppressWarnings("unchecked")
    public void testMatrixNormalization() throws NoSuchMethodException, InvocationTargetException, IllegalAccessException, NoSuchFieldException {
        testData[] data = new testData[]{
            new testData(
                new double[]{1,1,1,1}, new double[]{
                -0.842696629213483, 0.280898876404494, 0.280898876404494, 0.280898876404494, 0, 0,
                0.280898876404494, -1.12359550561798, 0.280898876404494, 0.280898876404494, 0.280898876404494, 0,
                0.280898876404494, 0.280898876404494, -1.12359550561798, 0.280898876404494, 0, 0.280898876404494,
                0.280898876404494, 0.280898876404494, 0.280898876404494, -0.842696629213483, 0, 0,
                0.280898876404494, 0.280898876404494, 0.280898876404494, 0.561797752808989, -1.40449438202247, 0,
                0.561797752808989, 0.280898876404494, 0.280898876404494, 0.280898876404494, 0, -1.40449438202247,
                }
            ),
            new testData(
                new double[]{4,3,2,1}, new double[]{
                -0.803921568627451, 0.133986928104575, 0.535947712418301, 0.133986928104575, 0, 0,
                0.133986928104575, -1.20588235294118, 0.133986928104575, 0.535947712418301, 0.401960784313726, 0,
                0.535947712418301, 0.133986928104575, -1.20588235294118, 0.133986928104575, 0, 0.401960784313726,
                0.133986928104575, 0.535947712418301, 0.133986928104575, -0.803921568627451, 0, 0,
                0.133986928104575, 0.26797385620915, 0.133986928104575, 0.669934640522876, -1.20588235294118, 0,
                0.669934640522876, 0.133986928104575, 0.26797385620915, 0.133986928104575, 0, -1.20588235294118,
                }
            )
        };
        for(testData test : data) {
            MethylationHKY substitutionModel = new MethylationHKY();
            substitutionModel.initByName(
                    "kappa", test.kappa,
                    "alpha", test.alpha,
                    "beta", test.beta,
                    "gamma", test.gamma
            );
            Method method = MethylationHKY.class.getDeclaredMethod("setupRateMatrix");
            method.setAccessible(true);
            method.invoke(substitutionModel);
            Field privateField = GeneralSubstitutionModel.class.getDeclaredField("rateMatrix");
            privateField.setAccessible(true);
            double[][] rateMatrix = (double[][]) privateField.get(substitutionModel);
            double[] flatMatrix = flattenMatrix(rateMatrix, 6, 6);
            //printMatrix(flatMatrix, 6, 6);
            Assert.assertArrayEquals(test.result, flatMatrix, 1e-10);
        }

    }


    private double[] flattenMatrix(double[][] matrix, int rows, int columns) {
        double[] flatMatrix = new double[rows * columns];
        for(int row=0; row < rows; row++){
            for(int col=0; col < columns; col++){
                flatMatrix[row * rows + col] = matrix[row][col];
            }
        }
        return flatMatrix;
    }
}
