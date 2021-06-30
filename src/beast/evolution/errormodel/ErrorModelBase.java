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
        double prob;
        int states = datatype.getStateCount();
        if (datatype.isAmbiguousCode(observedState)) {
            prob = getStatePartial(observedState, trueState);
        } else if (observedState == trueState) {
            prob =  1 - epsilon.getValue();
        } else {
            prob = epsilon.getValue() / (states - 1);
        }
        return prob;
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
