package beast.evolution.datatype;

import beast.core.Input;
import beast.core.parameter.RealParameter;

public class SiFitTernaryWithErrorSampled extends DataTypeWithErrorBase {

    int[][] x = {
            {0},
            {1},
            {2},
            {0, 1, 2},
    };

    public SiFitTernaryWithErrorSampled() {
        stateCount = 3;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "012" + MISSING_CHAR;
    }

    final public Input<RealParameter> alphaInput = new Input<>("alpha", "alpha parameter in SiFit Ternary model", Input.Validate.REQUIRED);
    final public Input<RealParameter> betaInput = new Input<>("beta", "beta parameter in SiFit Ternary model",  Input.Validate.REQUIRED);

    private RealParameter alpha;
    private RealParameter beta;

    protected double[][] errorMatrix;

    @Override
    public void initAndValidate() {
        // init base
        super.initAndValidate();
        // init error parameters
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
        if (isAmbiguousCode(observedState)) {
            return 1.0 / stateCount;
        } else {
            return errorMatrix[observedState][trueState];
        }
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
        return "sifitTernaryWithErrorSampled";
    }
}
