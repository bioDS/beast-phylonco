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
 * Truncated discrete Gaussian error model for copy number data.
 *
 * <p>Model assumptions:
 * <ul>
 *     <li>If true CN = 0, observed CN is deterministically 0</li>
 *     <li>If true CN = μ > 0, observed CN X follows a truncated discrete Gaussian
 *         over [0, nstate-1]: P(X = x | μ, σ) = exp(-(x-μ)²/(2σ²)) / Z(μ, σ, nstate)</li>
 *     <li>Z(μ, σ, nstate) is the normalization constant over x ∈ {0, 1, ..., nstate-1}</li>
 *     <li>σ is the measurement error standard deviation</li>
 *     <li>When σ → 0, all probability concentrates at the true value (no error).</li>
 * </ul>
 * <p> Truncation prevents sampling impossible copy numbers beyond nstate-1. </p>
 */
public class DiscreteGaussianErrorModel implements GenerativeDistribution<IntegerCharacterMatrix> {

    // Parameters
    Value<IntegerCharacterMatrix> alignment;
    Value<Double> sigma;  // Error standard deviation

    // Parameter names
    public final String alignmentParamName = "alignment";
    public final String sigmaParamName = "sigma";
    private int nstate;

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
        this.nstate = extractNstate(alignment);
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
                    "Applies a truncated discrete Gaussian error model to copy-number data. " +
                            "True CN = 0 is observed as 0; true CN = μ > 0 has " +
                            "P(obs=x) ∝ exp(-(x-μ)²/(2σ²)) truncated to [0, nstate-1]. " +
                            "σ controls measurement error; truncation prevents impossible values."
    )

    /**
     * Extract nstate from the CopyNumberBD model
     */
    private int extractNstate(Value<IntegerCharacterMatrix> alignment) {
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

        double sig = sigma.value();
        if (sig < 0) {
            throw new IllegalArgumentException("sigma must be non-negative.");
        }

        for (int i = 0; i < original.getTaxa().ntaxa(); i++) {
            String taxon = original.getTaxa().getTaxaNames()[i];

            for (int j = 0; j < original.nchar(); j++) {
                int trueState = original.getState(taxon, j);
                int observedState = sampleWithError(trueState, sig);
                observed.setState(i, j, observedState);
            }
        }

        return new RandomVariable<>("C", observed, this);
    }

    /**
     * Sampling rule for truncated discrete Gaussian:
     * <ul>
     *   <li>If true CN = 0 → return 0 deterministically (no measurement error at zero).</li>
     *   <li>If true CN = μ > 0 → sample from truncated discrete Gaussian over [0, nstate-1]
     *       centered at μ with SD σ, renormalized to account for truncation.</li>
     * </ul>
     *
     * @param trueCN The true copy number (μ)
     * @param sigma The measurement error standard deviation (σ)
     * @return The observed copy number after applying measurement error
     */
    private int sampleWithError(int trueCN, double sigma) {

        // True CN = 0 always observed as 0
        if (trueCN == 0) return 0;

        // If sigma very small, return true value (no error)
        if (sigma < 1e-6) return trueCN;

        // Compute probabilities for [0, nstate-1]
        double[] probs = new double[nstate];
        double sum = 0.0;

        for (int x = 0; x < nstate; x++) {
            double diff = x - trueCN;
            // Unnormalized: exp(-(x-μ)²/(2σ²))
            probs[x] = Math.exp(-(diff * diff) / (2.0 * sigma * sigma));
            sum += probs[x];
        }

        // Normalize
        for (int x = 0; x < nstate; x++) {
            probs[x] /= sum;
        }

        // Sample from discrete distribution
        double u = RandomUtils.getRandom().nextDouble();
        double cumulative = 0.0;
        for (int x = 0; x < nstate; x++) {
            cumulative += probs[x];
            if (u <= cumulative) {
                return x;
            }
        }

        // Fallback
        return Math.min(trueCN, nstate - 1);
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