package beast.evolution.datatype;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;

@Description("Binary error model from SiFit paper")
public class BinaryWithError extends Binary implements DataTypeWithError {

    final public Input<RealParameter> alphaInput = new Input<>("alpha", "alpha parameter in SiFit Binary model", Input.Validate.REQUIRED);
    final public Input<RealParameter> betaInput = new Input<>("beta", "beta parameter in SiFit Binary model",  Input.Validate.REQUIRED);

    private RealParameter alpha;
    private RealParameter beta;

    protected double[][] errorMatrix;

    public BinaryWithError() {
        super();
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();
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

    public double[] getProbabilities(int observedState) {
        double[] prob = new double[errorMatrix.length];
        for (int i = 0; i < errorMatrix.length; i++) {
            prob[i] = getProbability(observedState, i);
        }
        return prob;
    }

    @Override
    public String getTypeDescription() {
        return "binaryWithError";
    }
}
