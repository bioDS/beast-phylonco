package phylonco.lphy.evolution.copynumbermodel;

import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.simulator.RandomUtils;
import org.apache.commons.math3.distribution.PascalDistribution;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Map;
import java.util.TreeMap;

/**
 * Negative Binomial error model for copy number data.
 *
 * <p> This model assumes:
 * <ul>
 *     <li> If the true copy number is 0, the observed copy number is deterministically 0. </li>
 *     <li> If the true copy number is a > 0, the observed copy number X follows a Negative Binomial distribution: </li>
 *          X ~ NB(mean = a, variance = a + a²/k)
 * </ul>
 * <p> where k is the dispersion parameter.</p>
 * <p>
 * Parameterisation used:
 *      <li> r = a² / (variance - a) </li>
 *      <li> p = a / variance </li>
 * <p>
 * This corresponds to the Negative Binomial distribution:
 * <li> P(X = x) = Γ(x+r) / (x! Γ(r)) · (1-p)^x · p^r </li>
 */
public class NegativeBinomialErrorModel implements GenerativeDistribution<IntegerCharacterMatrix> {

    // Parameters
    Value<IntegerCharacterMatrix> alignment;
    Value<Double> dispersion;

    // Parameter names
    public final String alignmentParamName = "alignment";
    public final String dispersionParamName = "dispersion";

    RandomGenerator random;

    public NegativeBinomialErrorModel(
            @ParameterInfo(
                    name = alignmentParamName,
                    narrativeName = "true copy number alignment",
                    description = "The error-free (true) copy number alignment."
            )
            Value<IntegerCharacterMatrix> alignment,

            @ParameterInfo(
                    name = dispersionParamName,
                    narrativeName = "dispersion parameter",
                    description = "Negative binomial dispersion parameter k (Var = μ + μ²/k)."
            )
            Value<Double> dispersion
    ) {
        this.alignment = alignment;
        this.dispersion = dispersion;
        this.random = RandomUtils.getRandom();
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = new TreeMap<>();
        map.put(alignmentParamName, alignment);
        map.put(dispersionParamName, dispersion);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(alignmentParamName)) {
            alignment = value;
        } else if (paramName.equals(dispersionParamName)) {
            dispersion = value;
        } else {
            throw new RuntimeException("Unrecognized parameter name: " + paramName);
        }
    }

    @GeneratorInfo(
            name = "NegativeBinomialErrorModel",
            narrativeName = "negative binomial error model",
            verbClause = "has",
            category = GeneratorCategory.TAXA_ALIGNMENT,
            description =
                    "Applies a Negative Binomial error model to copy-number data. " +
                            "True CN = 0 is observed as 0; true CN > 0 is drawn from NB(mean = a, var = a + a²/k)."
    )

    @Override
    public RandomVariable<IntegerCharacterMatrix> sample() {

        IntegerCharacterMatrix original = alignment.value();
        IntegerCharacterMatrix observed = new IntegerCharacterMatrix(
                original.getTaxa(), original.nchar());

        double k = dispersion.value();
        if (k <= 0) {
            throw new IllegalArgumentException("Dispersion k must be positive.");
        }

        for (int i = 0; i < original.getTaxa().ntaxa(); i++) {
            String taxon = original.getTaxa().getTaxaNames()[i];

            for (int j = 0; j < original.nchar(); j++) {
                int trueState = original.getState(taxon, j);
                int observedState = sampleWithError(trueState, k);
                observed.setState(i, j, observedState);
            }
        }

        return new RandomVariable<>("C", observed, this);
    }

    /**
     * Sampling rule:
     * <ul>
     *   <li>If true CN = 0 → return 0.</li>
     *   <li>If true CN = a > 0 → sample from NB(mean = a, variance = a + a²/k).</li>
     * </ul>
     */
    private int sampleWithError(int trueCN, double k) {

        // True CN = 0 always observed as 0
        if (trueCN == 0) return 0;

        // Negative Binomial for true CN > 0
        double mean = (double) trueCN;
        double variance = mean + (mean * mean) / k;  // Var = μ + μ²/k

        // Convert to NB(r, p) parameterisation
        // mean = r(1-p)/p, var = mean + mean^2 / r
        double p = mean / variance;
        double r = (mean * mean) / (variance - mean);

        if (p <= 0.0 || p >= 1.0 || r <= 0.0) {
            // Fallback if parameters invalid
            return trueCN;
        }

        try {
            PascalDistribution pascalDist =
                    new PascalDistribution(random, (int) Math.round(r), p);
            int sample = pascalDist.sample();
            return Math.max(0, sample);

        } catch (Exception e) {
            // Fallback to true state if sampling fails
            return trueCN;
        }
    }

    public Value<IntegerCharacterMatrix> getAlignment() {
        return alignment;
    }

    public Value<Double> getDispersion() {
        return dispersion;
    }

    public IntegerCharacterMatrix getOriginalAlignment() {
        return alignment.value();
    }
}