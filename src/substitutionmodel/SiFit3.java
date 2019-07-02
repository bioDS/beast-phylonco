package substitutionmodel;

import beast.core.parameter.RealParameter;
import beast.core.Input;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;
import datatype.TernaryWithError;

import java.lang.reflect.InvocationTargetException;

public class SiFit3 extends GeneralSubstitutionModel {
    final public Input<RealParameter> lambdaDInput = new Input<>("lambdaD", "lambda deletions parameter in SiFit Ternary model",  Input.Validate.REQUIRED);
    final public Input<RealParameter> lambdaLInput = new Input<>("lambdaL", "lambda LOH parameter in SiFit Ternary model",  Input.Validate.REQUIRED);

    private RealParameter lambdaD;
    private RealParameter lambdaL;

    protected double[] frequencies;

    public SiFit3() {
        try {
            initAndValidate();
        } catch (Exception e) {
            e.printStackTrace();;
            throw new RuntimeException("initAndValidate() call failed when constructing SiFit3()");
        }
    }

    @Override
    public void initAndValidate() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
        frequenciesInput.setRule(Input.Validate.OPTIONAL);
        
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
        double[][] matrix = {
                {-1, 1, 0},
                {lambdaSum / 2, -lambdaSum, lambdaSum / 2},
                {0, lambdaD, -lambdaD}
        };
        rateMatrix = matrix;
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
        double[] freqs = {pi0, pi1, pi2};
        frequencies = freqs;
    }

    @Override
    protected double[][] getRateMatrix() {
        return rateMatrix.clone();
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
