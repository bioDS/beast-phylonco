package datatype;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.Binary;

// Binary error model from SiFit paper
public class BinaryWithError extends Binary implements DataTypeWithError {

    public Input<RealParameter> alphaInput = new Input<>("alpha", "alpha parameter in SiFit Binary model", Input.Validate.REQUIRED);
    public Input<RealParameter> betaInput = new Input<>("beta", "beta parameter in SiFit Binary model",  Input.Validate.REQUIRED);

    private RealParameter alpha;
    private RealParameter beta;

    protected double[][] errorMatrix;

    public BinaryWithError() {
        super();
        alpha = alphaInput.get();
        beta = betaInput.get();
        setupErrorMatrix();
    }

    @Override
    public void setupErrorMatrix() {
        setupErrorMatrix(alpha.getValue(), beta.getValue());
    }

    private void setupErrorMatrix(double alpha, double beta) {
        double[][] matrix = {
                {1 - alpha, beta},
                {alpha, 1 - beta}
        };
        errorMatrix = matrix;
    }

    public double getProbability(int observedState, int trueState) {
        return errorMatrix[observedState][trueState];
    }

    @Override
    public String getTypeDescription() {
        return "binaryWithError";
    }
}
