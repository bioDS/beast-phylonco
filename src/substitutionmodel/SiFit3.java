package substitutionmodel;

import beast.core.parameter.RealParameter;
import beast.core.Input;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.evolution.datatype.TernaryWithError;

import java.lang.reflect.InvocationTargetException;

/**
 * Implements the SiFit model of genotype substitution from Zafar et al. (2017)
 *
 * SiFit: inferring tumor trees from single-cell sequencing data under finite-sites models.
 * https://doi.org/10.1186/s13059-017-1311-2
 *
 * Q matrix
 *              0           1           2
 *	0 |        -1,          1,          0  |
 *	1 | (D+L) / 2,       -D-L,  (D+L) / 2  |
 *	2 |         0,          D,         -D  |
 *
 * with stationary distribution
 *
 *  x = 1 + 2 / (D+L) + 1/D
 *
 *  pi0 = 1 / (x)
 *  pi1 = 2 / (x * (D+L))
 *  pi2 = 1 / (x * D)
 */
public class SiFit3 extends GeneralSubstitutionModel {
    final public Input<RealParameter> lambdaDInput = new Input<>("lambdaD", "lambda deletions parameter in SiFit Ternary model",  Input.Validate.REQUIRED);
    final public Input<RealParameter> lambdaLInput = new Input<>("lambdaL", "lambda LOH parameter in SiFit Ternary model",  Input.Validate.REQUIRED);

    private RealParameter lambdaD;
    private RealParameter lambdaL;

    protected double[] frequencies;

    public SiFit3() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
        frequenciesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        if (ratesInput.get() != null) {
            throw new IllegalArgumentException("the rates attribute should not be used.");
        }
        if (frequenciesInput.get() != null) {
            throw new IllegalArgumentException("the frequencies should not be set, these can be calculated from the model parameters.");
        }

        lambdaD = lambdaDInput.get();
        lambdaL = lambdaLInput.get();
        updateMatrix = true;
        nrOfStates = 3;
        rateMatrix = new double[nrOfStates][nrOfStates];

        try {
            eigenSystem = createEigenSystem();
        } catch (ClassNotFoundException | InstantiationException | IllegalAccessException | InvocationTargetException e) {
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
                {-1, 1, 0},
                {lambdaSum / 2, -lambdaSum, lambdaSum / 2},
                {0, lambdaD, -lambdaD}
        };
    }

    protected void setupFrequencies() {
        setupFrequencies(lambdaD.getValue(), lambdaL.getValue());
    }

    private void setupFrequencies(double lambdaD, double lambdaL) {
        double lambdaSum = lambdaD + lambdaL;
        double x = 1 + 2 / (lambdaSum) + 1/lambdaD;
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
        return 3;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof TernaryWithError;
    }
}
