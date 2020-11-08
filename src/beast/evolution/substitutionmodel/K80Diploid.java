package beast.evolution.substitutionmodel;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.NucleotideDiploid;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;

public class K80Diploid extends GeneralSubstitutionModel {
    final public Input<RealParameter> kappaInput = new Input<>("kappa", "kappa the transition-transversion ratio in the Kimura model",  Input.Validate.REQUIRED);
    final public Input<RealParameter> lambdaLInput = new Input<>("lambdaL", "lambda L the combined rate of LOH and deletions in the SiFit model",  Input.Validate.REQUIRED);

    private RealParameter kappa;
    private RealParameter lambdaL;

    protected double[] frequencies;

    public K80Diploid() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
        frequenciesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        kappa = kappaInput.get();
        lambdaL = lambdaLInput.get();
        updateMatrix = true;
        nrOfStates = 10;
        rateMatrix = new double[nrOfStates][nrOfStates];

        try {
            eigenSystem = createEigenSystem();
        } catch(Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    protected void setupRelativeRates() {}

    @Override
    protected void setupRateMatrix() {
        setupFrequencies();
        setupRateMatrixUnnormalized();
        normalize();
    }

    @Override
    protected void setupRateMatrixUnnormalized() {
        setupRateMatrixUnnormalized(kappa.getValue(), lambdaL.getValue());
    }

    // instantaneous matrix Q
    private void setupRateMatrixUnnormalized(double kappa, double lambdaL) {
        rateMatrix = new double[][] {
                {0, 1, kappa, 1, 0, 0, 0, 0, 0, 0},
                {1, 0, 0.5, kappa / 2, 1, kappa / 2, 0.5, 0, 0, 0},
                {kappa, 0.5, 0, 0.5, 0, 0.5, 0, kappa, 0.5, 0},
                {1, kappa / 2, 0.5, 0, 0, 0, 0.5, 0, kappa / 2, 1},
                {0, 1, 0, 0, 0, 1, kappa, 0, 0, 0},
                {0, kappa / 2, 0.5, 0, 1, 0, 0.5, 1, kappa / 2, 0},
                {0, 0.5, 0, 0.5, kappa, 0.5, 0, 0, 0.5, kappa},
                {0, 0, kappa, 0, 0, 1, 0, 0, 1, 0},
                {0, 0, 0.5, kappa / 2, 0, kappa / 2, 0.5, 1, 0, 1},
                {0, 0, 0, 1, 0, 0, kappa, 0, 1, 0}
        };
        // deletion of one allele
        setupDeletion(rateMatrix, lambdaL);
        // calculate diagonal entries
        setupDiagonal(rateMatrix);
    }

    private void setupDeletion(double[][] rateMatrix, double lambdaL) {
        rateMatrix[1][0] += lambdaL; // AC to A-
        rateMatrix[2][0] += lambdaL; // AG to A-
        rateMatrix[3][0] += lambdaL; // AT to A-
        rateMatrix[1][4] += lambdaL; // AC to C-
        rateMatrix[5][4] += lambdaL; // CG to C-
        rateMatrix[6][4] += lambdaL; // CT to C-
        rateMatrix[2][7] += lambdaL; // AG to G-
        rateMatrix[5][7] += lambdaL; // CG to G-
        rateMatrix[8][7] += lambdaL; // GT to G-
        rateMatrix[3][9] += lambdaL; // AT to T-
        rateMatrix[6][9] += lambdaL; // CT to T-
        rateMatrix[8][9] += lambdaL; // GT to T-
    }

    private void setupDiagonal(double[][] rateMatrix) {
        for (int i = 0; i < nrOfStates; i++) {
            double sum = 0;
            for (int j = 0; j < nrOfStates; j++) {
                sum += rateMatrix[i][j];
            }
            rateMatrix[i][i] = -sum + rateMatrix[i][i];
        }
    }

    private void normalize() {
        double[] frequencies = getFrequencies();
        double f = 0.0;
        for (int i = 0; i < nrOfStates; i++) {
            f += frequencies[i] * -rateMatrix[i][i];
        }
        f = 1 / f;
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = f * rateMatrix[i][j];
            }
        }
    }

    protected void setupFrequencies() {
        double t = 10000;
        double[] matrix = new double[nrOfStates * nrOfStates];
        // get equilibrium frequencies using un-normalized matrix
        getTransitionProbabilities(null, t, 0, 1, matrix, false);
        frequencies = new double[nrOfStates];
        for (int i = 0; i < frequencies.length; i++) {
            frequencies[i] = matrix[i];
        }
    }

    @Override
    public double[] getFrequencies() {
        setupFrequencies();
        return frequencies;
    }

    @Override
    public int getStateCount() {
        return nrOfStates;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof NucleotideDiploid;
    }

}
