package beast.evolution.errormodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;

@Description("Error model base implementation")
public class ErrorModelBase extends ErrorModel {

    final public Input<RealParameter> epsilonInput = new Input<>("epsilon", "the per state error rate", Input.Validate.REQUIRED);

    private RealParameter epsilon;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        epsilon = epsilonInput.get();
    }

    @Override
    public double getProbability(int observedState, int trueState) {
        int states = datatype.getStateCount();
        if (datatype.isAmbiguousCode(observedState)) {
            return 1.0 / states;
        } else if (observedState == trueState) {
            return 1 - epsilon.getValue();
        } else {
            return epsilon.getValue() / (states - 1);
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
    public boolean canHandleDataType(DataType datatype) {
        return true;
    }
}
