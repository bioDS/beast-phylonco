package beast.evolution.substitutionmodel;

import beast.core.parameter.RealParameter;
import beast.core.Input;
import beast.evolution.datatype.Binary;
import beast.evolution.datatype.DataType;
import beast.evolution.substitutionmodel.GeneralSubstitutionModel;

/**
 * Implements a general binary substitution model
 */
public class BinarySubstitutionModel extends GeneralSubstitutionModel {
    final public Input<RealParameter> lambdaInput = new Input<>("lambda", "lambda the rate of deletion and back mutation",  Input.Validate.REQUIRED);

    private RealParameter lambda;

    protected double[] frequencies;

    public BinarySubstitutionModel() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
        frequenciesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        lambda = lambdaInput.get();
        updateMatrix = true;
        nrOfStates = 2;
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
        setupRateMatrix(lambda.getValue());
    }

    // instantaneous matrix Q
    private void setupRateMatrix(double lambda) {
        rateMatrix = new double[][] {
                {-1, 1},
                {lambda, -lambda}
        };
        normalize(rateMatrix);
    }

    private void normalize(double[][] rateMatrix) {
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
        setupFrequencies(lambda.getValue());
    }

    private void setupFrequencies(double lambda) {
        double pi0 = lambda / (lambda + 1);
        double pi1 = 1 / (lambda + 1);
        frequencies = new double[] {pi0, pi1};
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
        return dataType instanceof Binary;
    }
}
