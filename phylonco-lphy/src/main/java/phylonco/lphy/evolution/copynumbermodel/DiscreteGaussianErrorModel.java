package phylonco.lphy.evolution.copynumbermodel;

import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.simulator.RandomUtils;

import java.util.Map;
import java.util.TreeMap;

/**
 * Discrete Gaussian error model for copy number data.
 *
 * <p> This model assumes:
 * <ul>
 *     <li> If the true copy number is 0, the observed copy number is deterministically 0. </li>
 *     <li> If the true copy number is μ > 0, the observed copy number X follows a discrete Gaussian: </li>
 *          P(X = x | μ, σ) = exp(-(x-μ)²/(2σ²)) / Z(μ, σ)
 *     <li> where Z(μ, σ) is the normalization constant summed over all possible states </li>
 * </ul>
 * <p> σ is the measurement error standard deviation.</p>
 * <p> When σ → 0, all probability concentrates at the true value (no error). </p>
 * <p> The maximum copy number state is automatically inferred from the input alignment. </p>
 */
public class DiscreteGaussianErrorModel implements GenerativeDistribution<IntegerCharacterMatrix> {

    // Parameters
    Value<IntegerCharacterMatrix> alignment;
    Value<Double> sigma;  // Error standard deviation

    // Parameter names
    public final String alignmentParamName = "alignment";
    public final String sigmaParamName = "sigma";

    public DiscreteGaussianErrorModel(
            @ParameterInfo(
                    name = alignmentParamName,
                    narrativeName = "true copy number alignment",
                    description = "The error-free (true) copy number alignment."
            )
            Value<IntegerCharacterMatrix> alignment,

            @ParameterInfo(
                    name = sigmaParamName,
                    narrativeName = "measurement error SD",
                    description = "Standard deviation of Gaussian measurement error (σ). When σ=0, no error."
            )
            Value<Double> sigma
    ) {
        this.alignment = alignment;
        this.sigma = sigma;
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = new TreeMap<>();
        map.put(alignmentParamName, alignment);
        map.put(sigmaParamName, sigma);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(alignmentParamName)) {
            alignment = value;
        } else if (paramName.equals(sigmaParamName)) {
            sigma = value;
        } else {
            throw new RuntimeException("Unrecognized parameter name: " + paramName);
        }
    }

    @GeneratorInfo(
            name = "DiscreteGaussianErrorModel",
            narrativeName = "discrete Gaussian error model",
            verbClause = "has",
            category = GeneratorCategory.TAXA_ALIGNMENT,
            description =
                    "Applies a discrete Gaussian error model to copy-number data. " +
                            "True CN = 0 is observed as 0; true CN = μ > 0 has " +
                            "P(obs=x) ∝ exp(-(x-μ)²/(2σ²)). σ controls measurement error."
    )

    @Override
    public RandomVariable<IntegerCharacterMatrix> sample() {

        IntegerCharacterMatrix original = alignment.value();
        IntegerCharacterMatrix observed = new IntegerCharacterMatrix(
                original.getTaxa(), original.nchar());

        double sig = sigma.value();
        if (sig < 0) {
            throw new IllegalArgumentException("sigma must be non-negative.");
        }

        // Infer maximum copy number from the alignment data
        int maxCN = inferMaxCopyNumber(original);

        for (int i = 0; i < original.getTaxa().ntaxa(); i++) {
            String taxon = original.getTaxa().getTaxaNames()[i];

            for (int j = 0; j < original.nchar(); j++) {
                int trueState = original.getState(taxon, j);
                int observedState = sampleWithError(trueState, sig, maxCN);
                observed.setState(i, j, observedState);
            }
        }

        return new RandomVariable<>("C", observed, this);
    }

    /**
     * <p>Infer maximum copy number from the alignment.</p>
     * Finds the maximum state value across all taxa and sites,
     * then adds a buffer to allow for potential copy number gains during error sampling.
     */
    private int inferMaxCopyNumber(IntegerCharacterMatrix alignment) {
        int maxObserved = 0;

        for (int i = 0; i < alignment.getTaxa().ntaxa(); i++) {
            String taxon = alignment.getTaxa().getTaxaNames()[i];
            for (int j = 0; j < alignment.nchar(); j++) {
                int state = alignment.getState(taxon, j);
                if (state > maxObserved) {
                    maxObserved = state;
                }
            }
        }

        // Add buffer: allow observations up to maxObserved + 10 or at least 20
        // This ensures we have room for measurement error to push values higher
        return Math.max(maxObserved + 10, 20);
    }

    /**
     * Sampling rule:
     * <ul>
     *   <li>If true CN = 0 → return 0.</li>
     *   <li>If true CN = μ > 0 → sample from discrete Gaussian centered at μ with SD σ.</li>
     * </ul>
     */
    private int sampleWithError(int trueCN, double sigma, int maxCN) {

        // True CN = 0 always observed as 0
        if (trueCN == 0) return 0;

        // If sigma very small, return true value (no error)
        if (sigma < 1e-6) return trueCN;

        // Compute probabilities for all possible observed values
        double[] probs = new double[maxCN + 1];
        double sum = 0.0;

        for (int x = 0; x <= maxCN; x++) {
            double diff = x - trueCN;
            // Unnormalized: exp(-(x-μ)²/(2σ²))
            probs[x] = Math.exp(-(diff * diff) / (2.0 * sigma * sigma));
            sum += probs[x];
        }

        // Normalize
        for (int x = 0; x <= maxCN; x++) {
            probs[x] /= sum;
        }

        // Sample from discrete distribution
        double u = RandomUtils.getRandom().nextDouble();
        double cumulative = 0.0;
        for (int x = 0; x <= maxCN; x++) {
            cumulative += probs[x];
            if (u <= cumulative) {
                return x;
            }
        }

        // Fallback (shouldn't reach here)
        return trueCN;
    }

    public Value<IntegerCharacterMatrix> getAlignment() {
        return alignment;
    }

    public Value<Double> getSigma() {
        return sigma;
    }

    public IntegerCharacterMatrix getOriginalAlignment() {
        return alignment.value();
    }
}