package phylonco.lphy.evolution.copynumbermodel;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.simulator.RandomUtils;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.HashMap;
import java.util.Map;

/**
 * Copy Number Birth-Death Model for phylogenetic inference of copy number variation.
 *
 * <p>Models copy number evolution as a birth-death process where copy numbers can
 * gain (birth) or loss (death) by one unit per event.
 * State 0 is absorbing and cannot be recovered.</p>
 * <p>PARAMETERIZATION: Can be specified in two ways:
 * <ul>
 *   <li>Single rate: bdRate (assumes lambda=mu)</li>
 *   <li>Separate rates: lambda and mu (for future flexibility)</li>
 * </ul>
 * NOTE: Internally, lambda and mu are ALWAYS used for simulation.
 * When bdRate is provided, lambda and mu are set equal to bdRate.</p>
 *
 * <p>Used with {@link PhyloDiscrete} for simulating CNV data along phylogenetic trees.</p>
 */
public class CopyNumberBD extends DeterministicFunction<MarkovTraitEvolution<Integer>> implements MarkovTraitEvolution<Integer> {
    private Value<Double> lambda;
    private Value<Double> mu;
    private Value<Double> bdRate;
    private Value<Integer> rootState;
    private Value<Integer> nstate;

    public static final String birthRateString = "lambda";
    public static final String deathRateString = "mu";
    public static final String bdRateString = "bdRate";
    public static final String rootStateString = "rootState";
    public static final String nstateString = "nstate";

    // Constructor 1: Using single bdRate (lambda=mu)
    public CopyNumberBD(
            @ParameterInfo(name = bdRateString, description = "Birth-death rate when we assume lambda=mu")
            Value<Double> bdRate,
            @ParameterInfo(name = rootStateString, description = "Initial copy number at root", optional = true)
            Value<Integer> rootState,
            @ParameterInfo(name = nstateString, description = "Maximum number of copy number states (default=15)", optional = true)
            Value<Integer> nstate
    ) {
        super();
        assert (bdRate.value() >= 0) : "bdRate must be non-negative";
        this.rootState = rootState != null ? rootState : new Value<>("rootState", 2);
        this.nstate = nstate != null ? nstate : new Value<>("nstate", 15);
        this.bdRate = bdRate;
        this.lambda = new Value<>("lambda", bdRate.value());
        this.mu = new Value<>("mu", bdRate.value());
    }

    // Constructor 2: Using separate lambda and mu (allows lambda≠mu for future flexibility)
    public CopyNumberBD(
            @ParameterInfo(name = birthRateString, description = "Birth rate (gain)")
            Value<Double> lambda,
            @ParameterInfo(name = deathRateString, description = "Death rate (loss)")
            Value<Double> mu,
            @ParameterInfo(name = rootStateString, description = "Initial copy number at root", optional = true)
            Value<Integer> rootState,
            @ParameterInfo(name = nstateString, description = "Maximum number of copy number states (default=15)", optional = true)
            Value<Integer> nstate
    ) {
        super();
        // Use separate lambda and mu (allows lambda≠mu for future flexibility)
        assert (lambda.value() >= 0) : "lambda must be non-negative";
        assert (mu.value() >= 0) : "mu must be non-negative";
        this.rootState = rootState != null ? rootState : new Value<>("rootState", 2);
        this.nstate = nstate != null ? nstate : new Value<>("nstate", 15);
        this.lambda = lambda;
        this.mu = mu;
        this.bdRate = null;
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
        Map<String, Value> params = new HashMap<>();

        if (bdRate != null) {
            // User provided bdRate
            params.put(bdRateString, bdRate);
        } else {
            // User provided lambda and mu separately
            params.put(birthRateString, lambda);
            params.put(deathRateString, mu);
        }
        params.put(rootStateString, rootState);
        params.put(nstateString, nstate);
        return params;
    }

    public void setParam(String paramName, Value value) {
        if (paramName.equals(bdRateString)) {
            bdRate = value;
            // The simulator always uses lambda/mu internally
            lambda = new Value<>("lambda", bdRate.value());
            mu = new Value<>("mu", bdRate.value());
        } else if (paramName.equals(birthRateString)) {
            if (bdRate != null) {
                throw new RuntimeException("Cannot set lambda when using bdRate mode.");
            }
            lambda = value;
        } else if (paramName.equals(deathRateString)) {
            if (bdRate != null) {
                throw new RuntimeException("Cannot set mu when using bdRate mode.");
            }
            mu = value;
        } else if (paramName.equals(rootStateString)) {
            rootState = value;
        } else if (paramName.equals(nstateString)) {
            nstate = value;
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

    public Value<Double> getBdRate() {
        return bdRate;
    }
    public Value<Integer> getNstate() {
        return nstate;
    }
}