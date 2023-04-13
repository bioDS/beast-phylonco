package phylonco.beast.evolution.substitutionmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.datatype.Binary;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;

/**
 * Implements a general binary substitution model
 */
@Description("A binary substitution model with a single rate parameter")
public class BinarySubstitutionModel extends GeneralSubstitutionModel implements SubstitutionModel {
    final public Input<RealParameter> lambdaInput = new Input<>("lambda", "lambda the rate of deletion and back mutation",  Input.Validate.REQUIRED);

    private RealParameter lambda;

    protected double[] frequencies;

    public BinarySubstitutionModel() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
        frequenciesInput.setRule(Input.Validate.FORBIDDEN);
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
    public void setupRelativeRates() {}

    @Override
    public void setupRateMatrix() {
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
