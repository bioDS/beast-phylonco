package beast.evolution.errormodel;

import beast.core.CalculationNode;
import beast.core.Description;
import beast.core.Input;
import beast.evolution.datatype.DataType;

@Description("Error model abstract class")
public abstract class ErrorModel extends CalculationNode {

    final public Input<DataType> datatypeInput = new Input<>("datatype", "the datatype of the alignment, for example nucleotide etc", Input.Validate.REQUIRED);

    protected DataType datatype;

    protected double[][] errorMatrix;

    /**
     * returns the probability of observed state given the true state based on the error model
     * @param observedState index of observed state
     * @param trueState index of true state
     * @return probability of observed state given the true state
     */
    public abstract double getProbability(int observedState, int trueState);

    /**
     * returns a probability vector of the conditional probability of the observed state given each true state
     * @param observedState index of observed state
     * @return conditional probabilities of the observed state given the true state, for each possible true state
     */
    public abstract double[] getProbabilities(int observedState);

    /**
     * checks whether the error model can handle the input datatype
     * @param datatype the alignment datatype
     * @return true if the error model can handle the input datatype
     */
    public abstract boolean canHandleDataType(DataType datatype);

    /**
     * returns the state partial for ambiguous states, 1.0 if true state is in the set, or 0.0 otherwise
     * @param observedState observed state index
     * @param trueState true state index
     * @return state partial
     */
    public double getStatePartial(int observedState, int trueState) {
        double prob;
        boolean[] states = datatype.getStateSet(observedState);
        if (states[trueState]) {
            prob = 1.0;
        } else {
            prob = 0.0;
        }
        return prob;
    }

    @Override
    public void initAndValidate() {
        datatype = datatypeInput.get();
        if (datatype != null && !canHandleDataType(datatype)) {
            String message = "Error model cannot handle data type " + datatype.getTypeDescription();
            throw new IllegalArgumentException(message);
        }
    }

}
