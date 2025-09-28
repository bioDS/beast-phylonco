package phylonco.lphy.evolution.substitutionmodel;

import lphy.base.evolution.substitutionmodel.RateMatrix;
import lphy.core.model.Value;
import lphy.core.model.annotation.Citation;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import lphy.core.model.datatype.DoubleArray2DValue;

import java.util.Arrays;
import java.util.List;

/**
 * @author Alexei Drummond
 */
@Citation(
        value = "Alexey Kozlov, Joao Alves, Alexandros Stamatakis, David Posada (2021). CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data. bioRxiv 2020.07.31.230292.",
        title = "CellPhy: accurate and fast probabilistic inference of single-cell phylogenies from scDNA-seq data",
        authors = {"Kozlov", "Alves", "Stamatakis", "Posada"},
        year = 2021,
        DOI = "https://doi.org/10.1101/2020.07.31.230292"
)
public class GT10 extends RateMatrix {

    public static final String ratesParamName = "rates";
    public static final String freqParamName = "freq";

    public GT10(@ParameterInfo(name = ratesParamName, narrativeName = "relative rates", description = "the relative rates of the G10 process. Size 6.") Value<Double[]> rates,
                @ParameterInfo(name = freqParamName, narrativeName = "base frequencies", description = "the base frequencies of the G10 process. Size 10.") Value<Double[]> freq,
                @ParameterInfo(name = RateMatrix.meanRateParamName, narrativeName = "substitution rate", description = "the rate of substitution.", optional = true) Value<Number> meanRate) {

        super(meanRate);

        if (rates.value().length != 6) {
            throw new IllegalArgumentException("Rates must have 6 dimensions.");
        }

        if (freq.value().length != 10) {
            throw new IllegalArgumentException("Rates must have 10 dimensions.");
        }

        setParam(ratesParamName, rates);
        setParam(freqParamName, freq);
    }

    @GeneratorInfo(name = "gt10",
            verbClause = "is",
            narrativeName = "general time-reversible rate matrix on phased genotypes",
            description = "The GTR instantaneous rate matrix on unphased genotypes. Takes relative rates (6) and base frequencies (10) and produces an GT10 rate matrix.")
    public Value<Double[][]> apply() {
        Value<Double[]> rates = getRates();
        Value<Double[]> freq = getParams().get(freqParamName);
        return new DoubleArray2DValue(g10(rates.value(), freq.value()), this);
    }

    private Double[][] g10(Double[] rates, Double[] freqs) {

        int numStates = 10;

        Double[][] Q = new Double[numStates][numStates];

        String[] genotypes = {
                "AA", "AC", "AG", "AT", "CC",
                "CG", "CT", "GG", "GT", "TT"
        };


        // construct off-diagonals
        int rateIndex = 0;
        for (int i = 0; i < numStates; i++) {
            char fromParent1State = genotypes[i].charAt(0);
            char fromParent2State = genotypes[i].charAt(1);

            for (int j = i + 1; j < numStates; j++) {
                char toParent1State = genotypes[j].charAt(0);
                char toParent2State = genotypes[j].charAt(1);

                if (fromParent1State == toParent1State || fromParent2State == toParent2State || fromParent1State == toParent2State || fromParent2State == toParent1State) {

                    if (fromParent1State == toParent1State) {
                        rateIndex = getRateIndex(fromParent2State, toParent2State);
                    } else if (fromParent2State == toParent2State) {
                        rateIndex = getRateIndex(fromParent1State, toParent1State);
                    } else if (fromParent1State == toParent2State) {
                        rateIndex = getRateIndex(fromParent2State, toParent1State);
                    } else {
                        rateIndex = getRateIndex(fromParent1State, toParent2State);
                    }

                    Q[i][j] = rates[rateIndex] * freqs[j];
                    Q[j][i] = rates[rateIndex] * freqs[i];

                } else {
                    Q[i][j] = 0.0;
                    Q[j][i] = 0.0;
                }
            }
        }

        // construct diagonals
        for (
                int i = 0;
                i < numStates; i++) {
            double totalRate = 0.0;
            for (int j = 0; j < numStates; j++) {
                if (j != i) {
                    totalRate += Q[i][j];
                }
            }
            Q[i][i] = -totalRate;
        }

        // normalise rate matrix to one expected substitution per unit time
        normalize(freqs, Q, totalRateDefault1());

        return Q;
    }


    private int getRateIndex(char state1, char state2) {
        char[] states = { state1, state2 };
        Arrays.sort(states);
        String genotype = new String(states);
        return switch (genotype) {
            case "AC" -> 0; // r_AC
            case "AG" -> 1; // r_AG
            case "AT" -> 2; // r_AT
            case "CG" -> 3; // r_CG
            case "CT" -> 4; // r_CT
            case "GT" -> 5; // r_GT
            default -> throw new IllegalArgumentException("Invalid nucleotide pair: " + state1 + state2);
        };
    }




    public Value<Double[]> getRates() {
        return getParams().get(ratesParamName);
    }

    public Value<Double[]> getFreq() {
        return getParams().get(freqParamName);
    }
}