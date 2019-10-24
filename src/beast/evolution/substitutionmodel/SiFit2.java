package beast.evolution.substitutionmodel;

import beast.core.parameter.RealParameter;
import beast.core.Input;
import beast.evolution.datatype.Binary;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;

/**
 * Implements the two state SiFit model of genotype substitution from Zafar et al. (2017)
 *
 * SiFit: inferring tumor trees from single-cell sequencing data under finite-sites models.
 * https://doi.org/10.1186/s13059-017-1311-2
 *
 * Q matrix
 *              0           1
 *	0 |        -1,          1  |
 *	1 | (D+L) / 2, -(D+L) / 2  |
 *
 * with stationary distribution
 *
 *  x = 1 + (D+L) / 2
 *
 *  pi0 = (D+L) / (2 * x)
 *  pi1 = 1 / (x)
 */
public class SiFit2 extends GeneralSubstitutionModel {
    final public Input<RealParameter> lambdaDInput = new Input<>("lambdaD", "lambda D the rate of deletions in the SiFit Ternary model",  Input.Validate.REQUIRED);
    final public Input<RealParameter> lambdaLInput = new Input<>("lambdaL", "lambda L the rate of LOH in the SiFit Ternary model",  Input.Validate.REQUIRED);

    private RealParameter lambdaD;
    private RealParameter lambdaL;

    protected double[] frequencies;

    public SiFit2() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
        frequenciesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        lambdaD = lambdaDInput.get();
        lambdaL = lambdaLInput.get();
        updateMatrix = true;
        nrOfStates = 2;
        rateMatrix = new double[nrOfStates][nrOfStates];

        try {
            eigenSystem = createEigenSystem();
        } catch(Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    protected void setupRelativeRates() {}

    @Override
    protected void setupRateMatrix() {
        setupRateMatrix(lambdaD.getValue(), lambdaL.getValue());
    }

    // instantaneous matrix Q
    private void setupRateMatrix(double lambdaD, double lambdaL) {
        double lambdaSum = lambdaD + lambdaL;
        rateMatrix = new double[][] {
                {-1, 1},
                {lambdaSum / 2, -lambdaSum / 2}
        };
        normalize(rateMatrix);
    }

    private void normalize(double[][] rateMatrix) {
        double[] frequencies = getFrequencies();
        double f = 0.0;
        for (int i = 0; i < nrOfStates; i++) {
            f += frequencies[i] * -rateMatrix[i][i];
        }
        f = 1 / f;
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = f * rateMatrix[i][j];
            }
        }
    }

    protected void setupFrequencies() {
        setupFrequencies(lambdaD.getValue(), lambdaL.getValue());
    }

    private void setupFrequencies(double lambdaD, double lambdaL) {
        double lambdaSum = lambdaD + lambdaL;
        double x = 1 + lambdaSum / 2;
        double pi0 = lambdaSum / (2 * x);
        double pi1 = 1 / x;
        frequencies = new double[] {pi0, pi1};
    }

    @Override
    public double[] getFrequencies() {
        setupFrequencies();
        return frequencies;
    }

    @Override
    public int getStateCount() {
        return nrOfStates;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof Binary;
    }
}
