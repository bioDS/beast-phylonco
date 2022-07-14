package phylonco.beast.evolution.errormodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.Nucleotide;

/**
 * Implements the PMD error model from Rambaut et al. (2009)
 *
 * Accommodating the Effect of Ancient DNA Damage on Inferences of Demographic Histories
 * https://doi.org/10.1093/molbev/msn256
 *
 */
@Description("Error model base implementation")
public class PMDErrorModel extends ErrorModel {

    final public Input<RealParameter> decayInput = new Input<>("decayRate", "the rate of exponential decay rate r", Input.Validate.REQUIRED);

    private RealParameter decayRate;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        decayRate = decayInput.get();

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

    // add time t to superclass
    public double getProbability(int observedState, int trueState, double t) {
        double prob;
        int states = datatype.getStateCount();
        if (datatype.isAmbiguousCode(observedState)) {
            prob = getStatePartial(observedState, trueState);
        } else if (observedState == trueState) {
            prob =  Math.exp(-decayRate.getValue() * t);
        } else if isTransition(observedState, trueState) {
            prob = 1 - Math.exp(-decayRate.getValue() * t);
        } else {
            prob = 0
        }
        return prob;
    }

    private boolean isTransition(int i, int j) {
        String k = Nucleotide.getCharacter(i) + Nucleotide.getCharacter(j);
        if (k.equals("AG") || k.equals("GA")) {
            return true;
        } else if (k.equals("CT") || k.equals("TC")) {
            return true;
        } else{
            return false;
        }
    }

    // add time t to superclass
    public double[] getProbabilities(int observedState, double t) {
        int states = datatype.getStateCount();
        double[] prob = new double[states];
        for (int i = 0; i < states; i++) {
            prob[i] = getProbability(observedState, i, t);
        }
        return prob;
    }

    @Override
    public boolean canHandleDataType(DataType datatype) {
        return datatype instanceof Nucleotide;
    }
}