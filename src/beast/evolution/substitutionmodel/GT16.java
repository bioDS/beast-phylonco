package beast.evolution.substitutionmodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.datatype.DataType;
import beast.evolution.datatype.NucleotideDiploid16;

/**
 * Implements the GT16 model for diploid genotypes from Kozlov et al. (2021)
 *
 * CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data
 * https://doi.org/10.1101/2020.07.31.230292
 *
 */
@Description("GT16 diploid substitution model from CellPhy paper")
public class GT16 extends GeneralSubstitutionModel {

    final public Input<RealParameter> nucRatesInput = new Input<>("nucRates", "rate parameters for AC, AG, AT, CG, CT, GT");

    private RealParameter rates;

    public GT16() {
        super.ratesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        rates = nucRatesInput.get();
        rates.setInputValue("keys", "AC AG AT CG CT GT");

        if (super.ratesInput.get() != null) {
            throw new IllegalArgumentException("the rates attribute should not be used. Use nucRates instead.");
        }

        frequencies = frequenciesInput.get();
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
    protected void setupRelativeRates() { }

    @Override
    protected void setupRateMatrix() {
        setupFrequencies();
        setupRateMatrixUnnormalized();
        normalize();
    }

    @Override
    protected void setupRateMatrixUnnormalized() {
        double rateAC = rates.getValue("AC");
        double rateAG = rates.getValue("AG");
        double rateAT = rates.getValue("AT");
        double rateCG = rates.getValue("CG");
        double rateCT = rates.getValue("CT");
        double rateGT = rates.getValue("GT");

        setupRateMatrixUnnormalized(rateAC, rateAG, rateAT, rateCG, rateCT, rateGT);
    }

    // instantaneous matrix Q
    private void setupRateMatrixUnnormalized(double rateAC, double rateAG, double rateAT,
                                             double rateCG, double rateCT, double rateGT) {
        double[] pi = frequencies.getFreqs();
        rateMatrix = new double[nrOfStates][nrOfStates];
        int bases = 4;
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                double result = 0.0;
                int fromFirst = i / bases; // first allele in from state
                int fromSecond = i % bases; // second allele in from state
                int toFirst = j / bases; // first allele in to state
                int toSecond = j % bases; // second allele in to state
                if (i != j && (fromFirst == toFirst || fromSecond == toSecond)) {
                        int first, second;
                        if (fromFirst == toFirst) {
                            first = Math.min(fromSecond, toSecond);
                            second = Math.max(fromSecond, toSecond);
                        } else {
                            first = Math.min(fromFirst, toFirst);
                            second = Math.max(fromFirst, toFirst);
                        }
                        int orderedPair = first * 10 + second;
                        switch (orderedPair) {
                            case 1: // 01 - AC
                                result = rateAC; // A -> C or C -> A
                                break;
                            case 2: // 02 - AG
                                result = rateAG; // A -> G or G -> A
                                break;
                            case 3: // 03 - AT
                                result = rateAT; // A -> T or T -> A
                                break;
                            case 12: // 12 - CG
                                result = rateCG; // C -> G or G -> C
                                break;
                            case 13: // 13 - CT
                                result = rateCT; // C -> T or T -> C
                                break;
                            case 23: // 23 - GT
                                result = rateGT; // G -> T or T -> G
                                break;
                            default:
                                result = 0.0;
                                break;
                        }
                } else {
                    result = 0.0; // not reachable in single mutation or diagonal
                }
                rateMatrix[i][j] = result * pi[j];
            }
        }
        // calculate diagonal entries
        setupDiagonal(rateMatrix);
    }

    private void setupDiagonal(double[][] rateMatrix) {
        for (int i = 0; i < nrOfStates; i++) {
            double sum = 0;
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j)
                    sum += rateMatrix[i][j];
            }
            rateMatrix[i][i] = -sum;
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
    public int getStateCount() {
        return nrOfStates;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof NucleotideDiploid16;
    }

}
