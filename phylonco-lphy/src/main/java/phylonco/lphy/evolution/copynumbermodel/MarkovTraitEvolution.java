package phylonco.lphy.evolution.copynumbermodel;


/**
 * Interface for trait evolution following a Markov process.
 * @param <T> The trait type
 */
public interface MarkovTraitEvolution<T> {
    /**
     * Evolves a trait from a given state over a specified time period.
     * @param currentState The starting trait state
     * @param timeInterval The amount of evolutionary time
     * @return The evolved trait value
     */
    T evolveTraitOverTime(T currentState, double timeInterval);

    /**
     * Samples an ancestral (root) trait value from the root distribution.
     * @return The root trait value
     */
    T sampleAncestralTrait();

    /**
     * Calculates the probability of transitioning from one state to another.
     * @param startState The initial trait state
     * @param endState The final trait state
     * @param timeInterval The amount of evolutionary time
     * @return The probability of the transition
     */
    double transitionProbability(T startState, T endState, double timeInterval);
}