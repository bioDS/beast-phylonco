package beast.evolution.datatype;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;

@Description("Binary error model from SiFit paper with fixed error rates")
public class BinaryWithError extends Binary implements DataTypeWithError {

    final public Input<RealParameter> alphaInput = new Input<>("alpha", "false positive probability", Input.Validate.REQUIRED);
    final public Input<RealParameter> betaInput = new Input<>("beta", "false negative probability",  Input.Validate.REQUIRED);

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
        if (isAmbiguousCode(observedState)) {
            return 1.0 / stateCount;
        } else {
            return errorMatrix[observedState][trueState];
        }
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
