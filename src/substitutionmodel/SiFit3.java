package substitutionmodel;

import beast.core.parameter.RealParameter;
import beast.core.Input;
import beast.evolution.datatype.DataType;
import datatype.DataTypeWithError;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.evolution.tree.Node;

//import beast.core.CalculationNode;
//import beast.evolution.substitutionmodel.EigenDecomposition;
//import beast.evolution.substitutionmodel.SubstitutionModel;

public class SiFit3 extends GeneralSubstitutionModel {
    public Input<RealParameter> alphaInput = new Input<>("alpha", "alpha parameter in SiFit Ternary model", Input.Validate.REQUIRED);
    public Input<RealParameter> betaInput = new Input<>("beta", "beta parameter in SiFit Ternary model",  Input.Validate.REQUIRED);
    public Input<RealParameter> lambdaDInput = new Input<>("lambdaD", "lambda deletions parameter in SiFit Ternary model",  Input.Validate.REQUIRED);
    public Input<RealParameter> lambdaLInput = new Input<>("lambdaL", "lambda LOH parameter in SiFit Ternary model",  Input.Validate.REQUIRED);

    protected double[] rateMatrix;
    protected double alpha;
    protected double beta;
    protected double lambdaD;
    protected double lambdaL;

    @Override
    public void initAndValidate() {
        if (ratesInput.get() != null) {
            throw new IllegalArgumentException("the rates attribute should not be used.");
        }
        // set Frequency.estimateInput flag to false for uniform frequencies or disallow use of flag
        frequencies = frequenciesInput.get();
        updateMatrix = true;
        nrOfStates = 3;
        alpha = alphaInput.get().getValue();
        beta = betaInput.get().getValue();
        lambdaD = lambdaDInput.get().getValue();
        lambdaL = lambdaLInput.get().getValue();
        setupRateMatrix(lambdaD, lambdaL);
    }

    @Override
    protected void setupRateMatrix() {
        setupRateMatrix(lambdaD, lambdaL);
    }

    // instantaneous matrix Q
    private void setupRateMatrix(double lambdaD, double lambdaL) {
        double lambdaSum = lambdaD + lambdaL;
        double[] matrix = {
                -1, 1, 0,
                lambdaSum / 2, -lambdaSum, lambdaSum / 2,
                0, lambdaD, -lambdaD
        };
        rateMatrix = matrix;
    }

    @Override
    public double[] getRateMatrix(Node node) {
        return rateMatrix;
    }

    @Override
    public int getStateCount() {
        return 3;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof DataTypeWithError;
    }

//    @Override
//    protected void setupRelativeRates() {
//        double lambdaSum = lambdaD + lambdaL;
//        relativeRates[0] = 1;
//        relativeRates[1] = 0;
//        relativeRates[2] = lambdaSum / 2;
//        relativeRates[3] = lambdaSum / 2;
//        relativeRates[4] = 0;
//        relativeRates[5] = lambdaD;
//    }
}
