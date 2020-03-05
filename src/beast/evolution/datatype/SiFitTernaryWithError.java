package beast.evolution.datatype;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;

/**
 * Implements the SiFit error model for ternary data from Zafar et al. (2017)
 *
 * SiFit: inferring tumor trees from single-cell sequencing data under finite-sites models.
 * https://doi.org/10.1186/s13059-017-1311-2
 *
 * The three states
 *  0 = Homozygous reference     (2 non-mutant alleles)
 *  1 = Heterozygous             (1 mutant allele)
 *  2 = Homozygous non-reference (2 mutant alleles)
 *  ? = Missing data             (0, 1 or 2 mutant alleles)
 *
 * with error parameters
 *  A = False positive rate
 *  B = False negative rate
 *
 * and error matrix
 *                      0           1                 2
 *	0 |   1-A-(A * B / 2),      B / 2,                0  |
 *	1 |                 A,      1 - B,                0  |
 *	2 |        A * B / 2 ,      B / 2,                1  |
 *
 */
@Description("Ternary error model from SiFit paper")
public class SiFitTernaryWithError extends Ternary implements DataTypeWithError {
    final public Input<RealParameter> alphaInput = new Input<>("alpha", "alpha parameter in SiFit Ternary model", Input.Validate.REQUIRED);
    final public Input<RealParameter> betaInput = new Input<>("beta", "beta parameter in SiFit Ternary model",  Input.Validate.REQUIRED);

    private RealParameter alpha;
    private RealParameter beta;

    protected double[][] errorMatrix;

    public SiFitTernaryWithError() {
        super();
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        alpha = alphaInput.get();
        beta = betaInput.get();
        setupErrorMatrix();
    }

    @Override
    public void setupErrorMatrix() {
        setupErrorMatrix(alpha.getValue(), beta.getValue());
    }

    private void setupErrorMatrix(double alpha, double beta) {
        double[][] matrix = {
                {1 - alpha - (alpha * beta / 2), beta / 2, 0},
                {alpha, 1 - beta, 0},
                {alpha * beta / 2, beta / 2, 1}
        };
        errorMatrix = matrix;
    }

    public double getProbability(int observedState, int trueState) {
        // if observed state is ?
        // return new double[]{1.0 / 3, 1.0 / 3, 1.0 / 3};
        return errorMatrix[observedState][trueState];
    }

    public double[] getProbabilities(int observedState) {
        double[] prob = new double[errorMatrix.length];
        for (int i = 0; i < errorMatrix.length; i++) {
            prob[i] = getProbability(observedState, i);
        }
        return prob;
    }

    @Override
    public String getTypeDescription() {
        return "sifitTernaryWithError";
    }
}
