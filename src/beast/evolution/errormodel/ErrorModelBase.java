package beast.evolution.errormodel;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;

public class ErrorModelBase extends ErrorModel {

    final public Input<RealParameter> epsilonInput = new Input<>("epsilon", "the per state error rate", Input.Validate.REQUIRED);

    private RealParameter epsilon;

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        epsilon = epsilonInput.get();
        datatype = datatypeInput.get();
    }

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
        int states = datatype.getStateCount();
        double[] prob = new double[states];
        for (int i = 0; i < states; i++) {
            prob[i] = getProbability(observedState, i);
        }
        return prob;
    }

    @Override
    public void setupErrorMatrix() {
        int states = datatype.getStateCount();
        for (int i = 0; i < states; i++) {
            for (int j = 0; j < states; j++) {
                errorMatrix[i][j] = getProbability(i, j);
            }
        }
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return true;
    }
}
