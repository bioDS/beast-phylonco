package beast.evolution.datatype;

import beast.core.Input;
import beast.core.parameter.RealParameter;

public class NucleotideWithError extends Nucleotide implements DataTypeWithError{

    final public Input<RealParameter> epsilonInput = new Input<>("epsilon", "epsilon parameter for error", Input.Validate.REQUIRED);

    private RealParameter epsilon;

    protected double[][] errorMatrix;

    public NucleotideWithError() {
        super();
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        epsilon = epsilonInput.get();
        setupErrorMatrix();
    }

    @Override
    public void setupErrorMatrix() {
        setupErrorMatrix(epsilon.getValue());
    }

    private void setupErrorMatrix(double epsilon) {
        double[][] matrix = {
                {1 - 3 * epsilon, epsilon, epsilon, epsilon},
                {epsilon, 1 - 3 * epsilon, epsilon, epsilon},
                {epsilon, epsilon, 1 - 3 * epsilon, epsilon},
                {epsilon, epsilon, epsilon, 1 - 3 * epsilon}
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
        return "nucleotideWithError";
    }
}
