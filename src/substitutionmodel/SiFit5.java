package substitutionmodel;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Quinary;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;

public class SiFit5 extends GeneralSubstitutionModel {
    final public Input<RealParameter> lambdaDInput = new Input<>("lambdaD", "lambda D the rate of recurrent point mutation in SiFit model",  Input.Validate.REQUIRED);
    final public Input<RealParameter> lambdaLInput = new Input<>("lambdaL", "lambda L the combined rate of deletion and loss of heterozygosity in SiFit model",  Input.Validate.REQUIRED);

    private RealParameter lambdaD;
    private RealParameter lambdaL;

    protected double[] frequencies;

    public SiFit5() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
        frequenciesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        lambdaD = lambdaDInput.get();
        lambdaL = lambdaLInput.get();
        updateMatrix = true;
        nrOfStates = 5;
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
                {0, 0, 0, 0, 0},
                {lambdaL, -1 - lambdaL, 1, 0, 0},
                {lambdaL / 2, lambdaD / 2, -lambdaSum, lambdaL / 2, lambdaD / 2},
                {0, 0, 0, 0, 0},
                {0, 0, lambdaD, lambdaL, -lambdaSum}
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
        // equilibrium frequencies
        double pi0 = 0; // non-zero
        double pi1 = 0;
        double pi2 = 0;
        double pi3 = 1 - pi0; // non-zero
        double pi4 = 0;
        frequencies = new double[] {pi0, pi1, pi2, pi3, pi4};
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
        return dataType instanceof Quinary;
    }
}
