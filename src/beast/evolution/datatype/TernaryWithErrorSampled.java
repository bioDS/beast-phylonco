package beast.evolution.datatype;

import beast.core.Input;
import beast.core.parameter.RealParameter;

public class TernaryWithErrorSampled extends DataTypeWithErrorBase {
    final public Input<RealParameter> epsilonInput = new Input<>("epsilon", "epsilon parameter in ternary error model (total error probability)", Input.Validate.REQUIRED);

    private RealParameter epsilon;

    protected double[][] errorMatrix;

    public TernaryWithErrorSampled() {
        super();
    }

    @Override
    public void initAndValidate() {
        // init base
        super.initAndValidate();
        // init error parameters
        epsilon = epsilonInput.get();
        setupErrorMatrix();
    }

    @Override
    public void setupErrorMatrix() {
        setupErrorMatrix(epsilon.getValue());
    }

    private void setupErrorMatrix(double epsilon) {
        double[][] matrix = {
                {1 - epsilon, epsilon / 2, 0},
				{epsilon, 1 - epsilon, epsilon},
				{0, epsilon / 2, 1 - epsilon}
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
        return "ternaryWithErrorSampled";
    }
}
