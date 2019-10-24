package beast.evolution.datatype;

import beast.evolution.datatype.DataType;

public interface DataTypeWithError extends DataType {
    /**
     * returns the probability of observed state given the true state based on the error model
     * @param observedState index of observed state
     * @param trueState index of true state
     * @return probability of observed state given the true state
     */
    public double getProbability(int observedState, int trueState);

    public double[] getProbabilities(int observedState);

    /**
     * set up the error matrix according to the error model parameters
     */
    public void setupErrorMatrix();

}
