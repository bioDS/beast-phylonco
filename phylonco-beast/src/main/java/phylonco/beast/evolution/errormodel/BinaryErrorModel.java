package phylonco.beast.evolution.errormodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.datatype.Binary;
import beast.base.evolution.datatype.DataType;

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
                 errorMatrix[observedState][trueState] = getProbability(observedState, trueState);
             }
         }
    }

    @Override
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

    @Override
    public double[] getProbabilities(int observedState) {
        if (updateMatrix) {
            setupErrorMatrix();
            updateMatrix = false;
        }
        return errorMatrix[observedState];
    }

    @Override
    public boolean canHandleDataType(DataType datatype) {
        return datatype instanceof Binary;
    }
}
