package beast.evolution.errormodel;

import beast.core.parameter.RealParameter;

public class ErrorModelBase extends ErrorModel {

    private RealParameter epsilon;

    @Override
    public double getProbability(int observedState, int trueState) {
        if (observedState == trueState) {
            return 1 - epsilon.getValue();
        } else {
            return epsilon.getValue() / datatype.getStateCount();
        }
    }

    @Override
    public double[] getProbabilities(int observedState) {
        return new double[0];
    }

    @Override
    public void setupErrorMatrix() {
        double[][] matrix = {};
        errorMatrix = matrix;
    }
}
