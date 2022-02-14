package phylonco.beast.evolution.substitutionmodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Ternary;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;

/**
 * Implements the three state SiFit model of genotype substitution from Zafar et al. (2017)
 *
 * SiFit: inferring tumor trees from single-cell sequencing data under finite-sites models.
 * https://doi.org/10.1186/s13059-017-1311-2
 *
 * Q matrix
 *              0             1             2
 *	0 |        -1,            1,            0  |
 *	1 | (D+L) / 2,       -D - L,    (D+L) / 2  |
 *	2 |         0,            D,           -D  |
 *
 * with stationary distribution
 *
 *  x = 1 + 2 / (D+L) + (1/D)
 *
 *  pi0 = 1 / (x)
 *  pi1 = 2 / (x * (D+L))
 *  pi2 = 1 / (x * D)
 *
 * Modifications to the substitution model:
 * We include a modified version of the model that allows mutational pathways from state 1 to 0 and state 1 to 2
 * Setting 'mutationPathInput' to 'true' uses the following modifications
 *
 * Q matrix
 *                      0                     1                      2
 *	0 |                -1,                    1,                     0  |
 *	1 | (1/6) + (D+L) / 2,       -D - L - (3/2),     (1/2) + (D+L) / 2  |
 *	2 |                 0,                    D,                    -D  |
 *
 * with stationary distribution
 *
 *  x = (1/6) + (D+L) / 2
 *  y = (1/2) + (D+L) / 2
 *
 *  pi0 = x / (x + 1 + (y / D))
 *  pi1 = 1 / (x + 1 + (y / D))
 *  pi2 = y / (x + 1 + (y / D)) * D
 *
 */
@Description("SiFit ternary substitution model")
public class SiFit3 extends GeneralSubstitutionModel {
    final public Input<RealParameter> lambdaDInput = new Input<>("lambdaD", "lambda D the rate of deletions in the SiFit Ternary model",  Input.Validate.REQUIRED);
    final public Input<RealParameter> lambdaLInput = new Input<>("lambdaL", "lambda L the rate of LOH in the SiFit Ternary model",  Input.Validate.REQUIRED);
    final public Input<Boolean> mutationPathInput = new Input<>("mutationPath", "allow mutation paths from state 1 to 0 and state 1 to 2", false);

    private RealParameter lambdaD;
    private RealParameter lambdaL;
    private Boolean mutationPath;

    protected double[] frequencies;

    public SiFit3() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
        frequenciesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        lambdaD = lambdaDInput.get();
        lambdaL = lambdaLInput.get();
        mutationPath = mutationPathInput.get();
        updateMatrix = true;
        nrOfStates = 3;
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
        if (mutationPath) {
            setupRateMatrixWithMutation(lambdaD.getValue(), lambdaL.getValue());
        } else {
            setupRateMatrix(lambdaD.getValue(), lambdaL.getValue());
        }
    }

    private void setupRateMatrixWithMutation(double lambdaD, double lambdaL) {
        double lambdaSum = lambdaD + lambdaL;
        rateMatrix = new double[][] {
                {-1, 1, 0},
                {1.0 / 6 + lambdaSum / 2, -lambdaSum - 2.0 / 3, 0.5 + lambdaSum / 2},
                {0, lambdaD, -lambdaD}
        };
        normalize(rateMatrix);
    }

    // instantaneous matrix Q
    private void setupRateMatrix(double lambdaD, double lambdaL) {
        double lambdaSum = lambdaD + lambdaL;
        rateMatrix = new double[][] {
                {-1, 1, 0},
                {lambdaSum / 2, -lambdaSum, lambdaSum / 2},
                {0, lambdaD, -lambdaD}
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
        if (mutationPath) {
            setupFrequenciesWithMutation(lambdaD.getValue(), lambdaL.getValue());
        } else {
            setupFrequencies(lambdaD.getValue(), lambdaL.getValue());
        }
    }

    private void setupFrequenciesWithMutation(double lambdaD, double lambdaL) {
        double lambdaSum = lambdaD + lambdaL;
        double x = 1.0 / 6 + lambdaSum / 2;
        double y = 1.0 / 2 + lambdaSum / 2;
        double pi0 = x / (x + 1 + y / lambdaD);
        double pi1 = 1 / (x + 1 + y / lambdaD);
        double pi2 = y / (x + 1 + y / lambdaD) * lambdaD;
        frequencies = new double[] {pi0, pi1, pi2};
    }

    private void setupFrequencies(double lambdaD, double lambdaL) {
        double lambdaSum = lambdaD + lambdaL;
        double x = 1 + 2 / (lambdaSum) + 1 / lambdaD;
        double pi0 = 1 / x;
        double pi1 = 2 / (x * lambdaSum);
        double pi2 = 1 / (x * lambdaD);
        frequencies = new double[] {pi0, pi1, pi2};
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
        return dataType instanceof Ternary;
    }
}
