package beast.evolution.errormodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.Binary;
import beast.evolution.datatype.DataType;

@Description("Binary error model with parameters as false positive and false negative probabilities")
public class BinaryErrorModel extends ErrorModel {
	
    final public Input<RealParameter> alphaInput = new Input<>("alpha", "the false positive probability", Input.Validate.REQUIRED);
    final public Input<RealParameter> betaInput = new Input<>("beta", "the false negative probability",  Input.Validate.REQUIRED);

    private RealParameter alpha;
    private RealParameter beta;

    @Override
    public void initAndValidate() {
        // init base
        super.initAndValidate();
        // init error parameters
        alpha = alphaInput.get();
        beta = betaInput.get();
    }

     private void setupErrorMatrix(double alpha, double beta) {
        double[][] matrix = {
                {1 - alpha, beta},
                {alpha, 1 - beta}
        };
        errorMatrix = matrix;
    }

    public double getProbability(int observedState, int trueState) {
        double a = alpha.getValue();
        double b = beta.getValue();
        double prob;
        if (datatype.isAmbiguousCode(observedState)) {
            prob = getStatePartial(observedState, trueState);
        } else {
            int errorStates = observedState * 10 + trueState;
            switch (errorStates) {
                case 0: // 00
                    prob = 1 - a;
                    break;
                case 1: // 01
                    prob = b;
                    break;
                case 10: // 10
                    prob = a;
                    break;
                case 11: // 11
                    prob = 1 - b;
                    break;
                default:
                    prob = 0.0;
                    break;
            }
        }

        return prob;
    }

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
        return datatype instanceof Binary;
    }
}
