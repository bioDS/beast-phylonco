package phylonco.lphy.evolution.copynumbermodel;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.simulator.RandomUtils;
import org.apache.commons.math3.random.RandomGenerator;
import java.util.Map;


public class CopyNumberBD extends DeterministicFunction<MarkovTraitEvolution<Integer>> implements MarkovTraitEvolution<Integer> {
    private Value<Double> lambda;
    private Value<Double> mu;
    private Value<Integer> rootState;

    public static final String birthRateString = "lambda";
    public static final String deathRateString = "mu";
    public static final String rootStateString = "rootState";

    // Constructor with all parameters including optional rootState
    public CopyNumberBD(
            @ParameterInfo(name = birthRateString, description = "The birth rate >= 0 (generation) of copy number variants") Value<Double> lambda,
            @ParameterInfo(name = deathRateString, description = "The death rate >= 0 (loss) of copy number variants") Value<Double> mu,
            @ParameterInfo(name = rootStateString, description = "The initial copy number state at the root node") Value<Integer> rootState
    ) {
        super();
        this.lambda = lambda;
        this.mu = mu;
        this.rootState = rootState != null ? rootState : new Value<>("rootState", 2);

        // Birth rate(lambda) and Death rate(mu) are non-negative
        assert (lambda.value() >= 0);
        assert (mu.value() >= 0);
    }

    @Override
    public Integer evolveTraitOverTime(Integer currentState, double timeInterval) {
        // Simulate copy number evolution using birth-death process
        return simulateCopiesOnBranchBin(0.0, timeInterval, currentState);
    }

    @Override
    public Integer sampleAncestralTrait() {
        return rootState.value();
    }

    @Override
    public double transitionProbability(Integer startState, Integer endState, double timeInterval) {
        throw new UnsupportedOperationException("Analytical transition probabilities not implemented yet");
    }

    public int simulateCopiesOnBranchBin(double startTime, double endTime, int startCopies) {
        double tCurrent = startTime;
        int m = startCopies; // Current copy number

        // Simulate events until we reach the end time or copies go extinct
        while (tCurrent < endTime && m > 0) {
            // Calculate event rates
            double birthRate = lambda.value() * m;
            double deathRate = mu.value() * m;
            double totalRate = birthRate + deathRate;

            // If no events are possible, break
            if (totalRate == 0) {
                break;
            }

            // Sample time to next event from exponential distribution
            RandomGenerator randomGen = RandomUtils.getRandom();
            double tNext = -Math.log(randomGen.nextDouble()) / totalRate;

            // Check if we've reached branch end
            if (tCurrent + tNext > endTime) {
                break;
            }

            // Update current time
            tCurrent += tNext;

            // Determine event type
            double u = randomGen.nextDouble();
            if (u < birthRate / totalRate) {
                // Birth event (copy gain)
                m += 1;
        } else {
                // Death event (copy loss)
                m -= 1;
            }
        }
        return m;
    }


    @Override
    public Value<MarkovTraitEvolution<Integer>> apply() {
        return new Value<>(null, this, this);
    }

    public Map<String, Value> getParams() {
        return Map.of(
                birthRateString, lambda,
                deathRateString, mu,
                rootStateString, rootState
        );
    }

    public void setParam(String paramName, Value value) {
        if (paramName.equals(birthRateString)) {
            lambda = value;
        } else if (paramName.equals(deathRateString)) {
            mu = value;
        } else if (paramName.equals(rootStateString)) {
            rootState = value;
        } else {
            throw new RuntimeException("Unrecognised parameter name: " + paramName);
        }
    }

    public Value<Double> getLambda() {
        return lambda;
    }

    public Value<Double> getMu() {
        return mu;
    }

    public Value<Integer> getRootState() {
        return rootState;
    }

}
