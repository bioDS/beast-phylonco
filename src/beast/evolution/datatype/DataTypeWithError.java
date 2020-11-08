package beast.evolution.datatype;

import beast.core.Description;
import beast.evolution.datatype.DataType;

@Description("Data type and error model where the error parameters can be sampled or fixed")
public interface DataTypeWithError extends DataType {
    /**
     * returns the probability of observed state given the true state based on the error model
     * @param observedState index of observed state
     * @param trueState index of true state
     * @return probability of observed state given the true state
     */
    public double getProbability(int observedState, int trueState);

	/**
     * returns a probability vector of the conditional probability of the observed state given each true state
     * @param observedState index of observed state
     * @return conditional probabilities of the observed state given the true state, for each possible true state
     */
    public double[] getProbabilities(int observedState);

    /**
     * set up the error matrix according to the error model parameters
     */
    public void setupErrorMatrix();

}
