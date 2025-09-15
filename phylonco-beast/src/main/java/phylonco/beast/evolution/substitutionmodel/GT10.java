package phylonco.beast.evolution.substitutionmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.inference.parameter.RealParameter;
import phylonco.beast.evolution.datatype.NucleotideDiploid10;

import java.util.Arrays;
import java.util.List;

/**
 * Implements the GT10 model for diploid genotypes from Kozlov et al. (2021)
 *
 * CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data
 * https://doi.org/10.1101/2020.07.31.230292
 *
 */
@Description("GT10 diploid substitution model from CellPhy paper")
public class GT10 extends GeneralSubstitutionModel implements SubstitutionModel {

    final public Input<RealParameter> nucRatesInput = new Input<>("nucRates", "rate parameters for AC, AG, AT, CG, CT, GT");

    private RealParameter rates;

    public GT10() {
        super.ratesInput.setRule(Input.Validate.OPTIONAL);
    }

    @Override
    public void initAndValidate() {
        rates = nucRatesInput.get();

        // validation checks
        if (super.ratesInput.get() != null) {
            throw new IllegalArgumentException("The rates attribute should not be used, use nucRates instead.");
        }

        if (rates == null) {
            throw new IllegalArgumentException("nucRates attribute needs to be specified.");
        } else if (rates.getDimension() != 6) {
            throw new IllegalArgumentException("nucRates dimension not equal to 6.");
        } else {
            List<String> keys = Arrays.asList("AC", "AG", "AT", "CG", "CT", "GT");
            for (String k: keys) {
                if (rates.getValue(k) == null)
                    throw new IllegalArgumentException("nucRates key needs to be specified for " + k);
            }
        }

        frequencies = frequenciesInput.get();
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
    public void setupRelativeRates() { }

    @Override
    public void setupRateMatrix() {
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
        String[] genotypes = {
                "AA", "AC", "AG", "AT", "CC",
                "CG", "CT", "GG", "GT", "TT"
        };
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                double result = 0.0;
                int fromFirst = baseIndex(genotypes[i].charAt(0)); // first allele in from state
                int fromSecond = baseIndex(genotypes[i].charAt(1)); // second allele in from state
                int toFirst = baseIndex(genotypes[j].charAt(0)); // first allele in to state
                int toSecond = baseIndex(genotypes[j].charAt(1)); // second allele in to state
                if (i != j && (fromFirst == toFirst || fromSecond == toSecond || fromFirst == toSecond || fromSecond == toFirst)) {
                    int first, second;
                    if (fromFirst == toFirst) {
                        first = Math.min(fromSecond, toSecond);
                        second = Math.max(fromSecond, toSecond);
                    } else if (fromSecond == toSecond) {
                        first = Math.min(fromFirst, toFirst);
                        second = Math.max(fromFirst, toFirst);
                    } else if (fromFirst == toSecond) {
                        first = Math.min(fromSecond, toFirst);
                        second = Math.max(fromSecond, toFirst);
                    } else {
                        first = Math.min(fromFirst, toSecond);
                        second = Math.max(fromFirst, toSecond);
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

    private int baseIndex(char base) {
        return switch (base) {
            case 'A' -> 0;
            case 'C' -> 1;
            case 'G' -> 2;
            case 'T' -> 3;
            default -> throw new IllegalArgumentException("Unknown base: " + base);
        };
    }


    @Override
    public int getStateCount() {
        return nrOfStates;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof NucleotideDiploid10;
    }

}
