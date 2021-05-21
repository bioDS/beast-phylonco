package beast.evolution.substitutionmodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.NucleotideDiploid16;

@Description("GT16 diploid substitution model from CellPhy paper")
public class GT16 extends GeneralSubstitutionModel {
    final public Input<RealParameter> rateACInput = new Input<>("rateAC", "the rate of A to C",  Input.Validate.REQUIRED);
    final public Input<RealParameter> rateAGInput = new Input<>("rateAG", "the rate of A to G",  Input.Validate.REQUIRED);
    final public Input<RealParameter> rateATInput = new Input<>("rateAT", "the rate of A to T",  Input.Validate.REQUIRED);
    final public Input<RealParameter> rateCGInput = new Input<>("rateCG", "the rate of C to G",  Input.Validate.REQUIRED);
    final public Input<RealParameter> rateCTInput = new Input<>("rateCT", "the rate of C to T",  Input.Validate.REQUIRED);
    final public Input<RealParameter> rateGTInput = new Input<>("rateGT", "the rate of G to T",  Input.Validate.REQUIRED);


    private RealParameter rateAC; // rate AC
    private RealParameter rateAG; // rate AG
    private RealParameter rateAT; // rate AT
    private RealParameter rateCG; // rate CG
    private RealParameter rateCT; // rate CT
    private RealParameter rateGT; // rate GT

    public GT16() {
        ratesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        if (ratesInput.get() != null) {
            throw new IllegalArgumentException("the rates attribute should not be used. Use the individual rates rateAC, rateCG, etc, instead.");
        }

        frequencies = frequenciesInput.get();

        rateAC = rateACInput.get();
        rateAG = rateAGInput.get();
        rateAT = rateATInput.get();
        rateCG = rateCGInput.get();
        rateCT = rateCTInput.get();
        rateGT = rateGTInput.get();

        updateMatrix = true;
        nrOfStates = 16;
        rateMatrix = new double[nrOfStates][nrOfStates];

        try {
            eigenSystem = createEigenSystem();
        } catch(Exception e) {
            e.printStackTrace();
        }
    }

    @Override
    protected void setupRelativeRates() {
        relativeRates[0] = rateAC.getValue(); // AA -> AC
        relativeRates[1] = rateAG.getValue(); // AA -> AG
        relativeRates[2] = rateAT.getValue(); // AA -> AT
    }

    @Override
    protected void setupRateMatrix() {
        setupFrequencies();
        setupRateMatrixUnnormalized();
        normalize();
    }

    @Override
    protected void setupRateMatrixUnnormalized() {
        setupRateMatrixUnnormalized(rateAC.getValue(), rateAG.getValue(), rateAT.getValue(),
                rateCG.getValue(), rateCT.getValue(), rateGT.getValue());
    }

    // instantaneous matrix Q
    private void setupRateMatrixUnnormalized(double rateAC, double rateAG, double rateAT,
                                             double rateCG, double rateCT, double rateGT) {
        double[] pi = frequencies.getFreqs();
        rateMatrix = new double[][] {

        };
        // calculate diagonal entries
        setupDiagonal(rateMatrix);
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
        frequencies.initAndValidate();
    }

    @Override
    public double[] getFrequencies() {
        return frequencies.getFreqs();
    }

    @Override
    public int getStateCount() {
        return nrOfStates;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof NucleotideDiploid16;
    }

}
