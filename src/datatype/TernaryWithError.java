package datatype;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;

// Ternary error model from SiFit paper
public class TernaryWithError extends DataType.Base implements DataTypeWithError {
    final public Input<RealParameter> alphaInput = new Input<>("alpha", "alpha parameter in SiFit Ternary model", Input.Validate.REQUIRED);
    final public Input<RealParameter> betaInput = new Input<>("beta", "beta parameter in SiFit Ternary model",  Input.Validate.REQUIRED);

    private RealParameter alpha;
    private RealParameter beta;

    protected double[][] errorMatrix;

    int[][] x = {
            {0},  // 0 Homozygous reference (wildtype)
            {1},  // 1 Heterozygous
            {2},  // 2 Homozygous non reference
            {0, 1, 2}, // ?
    };

    public TernaryWithError() {
        stateCount = 3;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "012" + MISSING_CHAR;

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
                {1 - alpha - (alpha * beta / 2), alpha, alpha * beta / 2},
                {beta / 2, 1 - beta, beta / 2},
                {0, 0, 1}
        };
        errorMatrix = matrix;
    }

    public double getProbability(int observedState, int trueState) {
        return errorMatrix[observedState][trueState];
    }

    @Override
    public String getTypeDescription() {
        return "ternaryWithError";
    }
}
