package phylonco.beast.evolution.errormodel;

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

        if (updateMatrix) {
            setupErrorMatrix();
            updateMatrix = false;
        }
    }

    @Override
    public void setupErrorMatrix() {
        if (errorMatrix == null) {
            errorMatrix = new double[datatype.mapCodeToStateSet.length][datatype.getStateCount()];
        }
        for (int trueState = 0; trueState < datatype.getStateCount(); trueState++) {
            for (int observedState = 0; observedState < datatype.mapCodeToStateSet.length; observedState++) {
                // rows are observed states X, columns are true states Y
                errorMatrix[observedState][trueState] = getProbability(observedState, trueState, 0.0);
            }
        }
    }

    @Override
    public double getProbability(int observedState, int trueState, double t) {
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
    public double[] getProbabilities(int observedState, double t) {
        int states = datatype.getStateCount();
        double[] prob = new double[states];
        for (int i = 0; i < states; i++) {
            prob[i] = getProbability(observedState, i, 0.0);
        }
        return prob;
    }

    @Override
    public boolean canHandleDataType(DataType datatype) {
        return true;
    }
}
