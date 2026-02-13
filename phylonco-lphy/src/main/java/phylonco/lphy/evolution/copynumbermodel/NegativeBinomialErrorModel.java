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
 * Truncated Negative Binomial error model for copy number data.
 *
 * <p>Model assumptions:
 * <ul>
 *     <li>If true CN = 0, observed CN is deterministically 0</li>
 *     <li>If true CN = a > 0, observed CN ~ NB(mean = a, variance = a + a²/k)
 *         truncated to [0, nstate-1], where k is the dispersion parameter</li>
 * </ul>
 *
 * <p>Parameterization: r = k, p = k/(a+k), giving
 * P(X = x) = Γ(x+r) / (x! Γ(r)) · (1-p)^x · p^r
 *
 * <p><strong>IMPORTANT:</strong> This model is TRUNCATED to the state space [0, nstate-1].
 * <ul>
 *     <li>Eg: nstate = 20 means maximum CN = 19 (states: 0,1,2,...,19)</li>
 * </ul>
 * All probability mass is renormalized to prevent
 * sampling impossible copy numbers.
 */
public class NegativeBinomialErrorModel implements GenerativeDistribution<IntegerCharacterMatrix> {

    // Parameters
    Value<IntegerCharacterMatrix> alignment;
    Value<Double> dispersion;

    // Parameter names
    public final String alignmentParamName = "alignment";
    public final String dispersionParamName = "dispersion";
    private int nstate;
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
        this.nstate = extractNstate(alignment);
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
    /**
     * Extract nstate from the CopyNumberBD model used to generate the alignment
     */
    private int extractNstate(Value<IntegerCharacterMatrix> alignment) {
        // Get the PhyloDiscrete generator
        if (alignment.getGenerator() instanceof PhyloDiscrete) {
            PhyloDiscrete phyloDiscrete = (PhyloDiscrete) alignment.getGenerator();
            Value<?> modelValue = phyloDiscrete.getModel();

            if (modelValue.value() instanceof CopyNumberBD) {
                CopyNumberBD cnvModel = (CopyNumberBD) modelValue.value();
                return cnvModel.getNstate().value();
            }
        }
        throw new IllegalStateException("Cannot extract nstate from alignment generator");
    }

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
     * Sampling rule for truncated Negative Binomial:
     * <ul>
     *   <li>If true CN = 0 → return 0 deterministically (no measurement error at zero).</li>
     *   <li>If true CN = a > 0 → sample from truncated NB(mean = a, variance = a + a²/k)
     *       over [0, nstate-1], renormalized to account for truncation.</li>
     * </ul>
     *
     * @param trueCN The true copy number (a)
     * @param k The dispersion parameter (controls overdispersion)
     * @return The observed copy number after applying measurement error
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
            // Build truncated distribution: compute probs for [0, nstate-1]
            double[] probs = new double[nstate];
            double sum = 0.0;

            PascalDistribution pascalDist =
                    new PascalDistribution(random, (int) Math.round(r), p);

            for (int x = 0; x < nstate; x++) {
                probs[x] = pascalDist.probability(x);
                sum += probs[x];
            }

            // Normalize
            for (int x = 0; x < nstate; x++) {
                probs[x] /= sum;
            }

            // Sample from truncated distribution
            double u = random.nextDouble();
            double cumulative = 0.0;
            for (int x = 0; x < nstate; x++) {
                cumulative += probs[x];
                if (u <= cumulative) {
                    return x;
                }
            }

            // Fallback
            return Math.min(trueCN, nstate - 1);

        } catch (Exception e) {
            return Math.min(trueCN, nstate - 1);
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