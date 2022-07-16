package phylonco.beast.evolution.errormodel;

import beast.evolution.datatype.DataType;
import beast.evolution.datatype.NucleotideDiploid16;
import org.junit.Test;

import java.util.List;

import static junit.framework.Assert.assertEquals;

public class GT16ErrorModelTest {

    private static double DELTA = 1e-10;

    /***
     * Given the true state is y, the conditional probability over all possible states is 1.
     * sum( P(s | y) ) = 1,  where s are all possible states.
     * Tests sum( P(s | y) ) = 1 for all possible true states.
     ***/
    @Test
    public void testGT16ErrorSumsToOne() {
        NucleotideDiploid16 datatype = new NucleotideDiploid16();

        GT16ErrorModel errorModel = new GT16ErrorModel();
        errorModel.initByName(
                "epsilon", "0.1",
                "delta", "0.2",
                "datatype", datatype
        );
        errorModel.initAndValidate();

        for (int trueState = 0; trueState < datatype.getStateCount(); trueState++) {
            double sum = 0.0;
            for (int observedState = 0; observedState < datatype.getStateCount(); observedState++) {
                sum += errorModel.getProbability(observedState, trueState, 0.0);
            }
            assertEquals(1.0, sum, DELTA);
        }

    }

    public double[][] getExpectedMatrix(double epsilon, double delta) {
        // Y is homozygous
        double a = 1 - epsilon + (1/2.0) * delta * epsilon;
        double b = (1 - delta) * (1/6.0) * epsilon;
        double c = (1/6.0) * delta * epsilon;
        // Y is heterozygous
        double d = (1/2.0) * delta + (1/6.0) * epsilon - (1/3.0) * delta * epsilon;
        double e = (1/6.0) * delta * epsilon;
        double f = (1 - delta) * (1/6.0) * epsilon;
        double g = (1 - delta) * (1 - epsilon);
        // rows are true states Y, columns are observed states X
        // the entries in the matrix are P(X | Y)
        double[][] expectedMatrix = {
                //      AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT
                {a, b, b, b, b, c, 0, 0, b, 0, c, 0, b, 0, 0, c}, // AA
                {d, g, f, f, 0, d, 0, 0, 0, f, e, 0, 0, f, 0, e}, // AC
                {d, f, g, f, 0, e, f, 0, 0, 0, d, 0, 0, 0, f, e}, // AG
                {d, f, f, g, 0, e, 0, f, 0, 0, e, f, 0, 0, 0, d}, // AT
                {d, 0, 0, 0, g, d, f, f, f, 0, e, 0, f, 0, 0, e}, // CA
                {c, b, 0, 0, b, a, b, b, 0, b, c, 0, 0, b, 0, c}, // CC
                {e, 0, f, 0, f, d, g, f, 0, 0, d, 0, 0, 0, f, e}, // CG
                {e, 0, 0, f, f, d, f, g, 0, 0, e, f, 0, 0, 0, d}, // CT
                {d, 0, 0, 0, f, e, 0, 0, g, f, d, f, f, 0, 0, e}, // GA
                {e, f, 0, 0, 0, d, 0, 0, f, g, d, f, 0, f, 0, e}, // GC
                {c, 0, b, 0, 0, c, b, 0, b, b, a, b, 0, 0, b, c}, // GG
                {e, 0, 0, f, 0, e, 0, f, f, f, d, g, 0, 0, 0, d}, // GT
                {d, 0, 0, 0, f, e, 0, 0, f, 0, e, 0, g, f, f, d}, // TA
                {e, f, 0, 0, 0, d, 0, 0, 0, f, e, 0, f, g, f, d}, // TC
                {e, 0, f, 0, 0, e, f, 0, 0, 0, d, 0, f, f, g, d}, // TG
                {c, 0, 0, b, 0, c, 0, b, 0, 0, c, b, b, b, b, a}  // TT
        };

        return expectedMatrix;
    }

    @Test
    public void testGT16SmallErrorProbabilities() {
        double epsilon = 0.01;
        double delta = 0.02;

        double [][] expectedMatrix = getExpectedMatrix(epsilon, delta);

        NucleotideDiploid16 datatype = new NucleotideDiploid16();

        GT16ErrorModel errorModel = new GT16ErrorModel();
        errorModel.initByName(
                "epsilon", Double.toString(epsilon),
                "delta", Double.toString(delta),
                "datatype", datatype
        );
        errorModel.initAndValidate();

        for (int trueState = 0; trueState < datatype.getStateCount(); trueState++) {
            for (int observedState = 0; observedState < datatype.getStateCount(); observedState++) {
                double calcProb = errorModel.getProbability(observedState, trueState, 0.0);
                double expectedProb = expectedMatrix[trueState][observedState];
                assertEquals(expectedProb, calcProb, DELTA);
            }
        }
    }

    @Test
    public void testGT16LargeErrorProbabilities() {
        double epsilon = 0.1;
        double delta = 0.5;

        double [][] expectedMatrix = getExpectedMatrix(epsilon, delta);

        NucleotideDiploid16 datatype = new NucleotideDiploid16();

        GT16ErrorModel errorModel = new GT16ErrorModel();
        errorModel.initByName(
                "epsilon", Double.toString(epsilon),
                "delta", Double.toString(delta),
                "datatype", datatype
        );
        errorModel.initAndValidate();

        for (int trueState = 0; trueState < datatype.getStateCount(); trueState++) {
            for (int observedState = 0; observedState < datatype.getStateCount(); observedState++) {
                double calcProb = errorModel.getProbability(observedState, trueState, 0.0);
                double expectedProb = expectedMatrix[trueState][observedState];
                assertEquals(expectedProb, calcProb, DELTA);
            }
        }
    }

    /***
     * Tests P(? | y) = 1 for all possible true states y.
     * Tests P(- | y) = 1 for all possible true states y.
     */

    @Test
    public void testGT16GapMissingPartialsEqualOne() {
        NucleotideDiploid16 datatype = new NucleotideDiploid16();

        GT16ErrorModel errorModel = new GT16ErrorModel();
        errorModel.initByName(
                "epsilon", "0.1",
                "delta", "0.2",
                "datatype", datatype
        );
        errorModel.initAndValidate();

        String ambiguousStr = "" + DataType.GAP_CHAR + DataType.MISSING_CHAR;
        List<Integer> codes = datatype.stringToEncoding(ambiguousStr);
        for (int observedState: codes) {
            for (int trueState = 0; trueState < datatype.getStateCount(); trueState++) {
                double expected = 1.0;
                double delta = 0.0;
                double modelProb = errorModel.getProbability(observedState, trueState, 0.0);
                assertEquals(expected, modelProb, delta);
            }
        }
    }

    /***
     * Given the true state is y, the conditional probability over all possible states is 1.
     * sum( P(s | y) ) = 1, where s are all possible states.
     * An ambiguous code with two states z = {s1, s2},
     * has probability P(z | t) = P(s1 | t) + P(s2 | t).
     * P(z | y) + sum( P(x | y) ) = 1, where x are all possible states excluding z.
     * Tests that P(z | y) + sum( P(x | y) ) = 1 for all possible true states.
     */
    @Test
    public void testGT16AmbiguousPartialsSumsToOne() {
        NucleotideDiploid16 datatype = new NucleotideDiploid16();

        GT16ErrorModel errorModel = new GT16ErrorModel();
        errorModel.initByName(
                "epsilon", "0.1",
                "delta", "0.2",
                "datatype", datatype
        );
        errorModel.initAndValidate();

        String ambiguousStr = "" + DataType.GAP_CHAR + DataType.MISSING_CHAR;
        for (int observedState = 0; observedState < datatype.mapCodeToStateSet.length; observedState++) {
            int[] observedArr = {observedState};
            String observedString = datatype.encodingToString(observedArr);
            if (datatype.isAmbiguousCode(observedState) && !ambiguousStr.contains(observedString)) {
                for (int trueState = 0; trueState < datatype.getStateCount(); trueState++) {
                    double sum = 0.0;
                    int[] observedCodes = datatype.getStatesForCode(observedState);
                    for (int otherStates = 0; otherStates < datatype.getStateCount(); otherStates++) {
                        // assumes only two ambiguous state codes
                        if (otherStates != observedCodes[0] && otherStates != observedCodes[1]) {
                            double modelProb = errorModel.getProbability(otherStates, trueState, 0.0);
                            sum += modelProb;
                        }
                    }
                    double expected = 1.0;
                    double modelProb = errorModel.getProbability(observedState, trueState, 0.0);
                    sum += modelProb;
                    assertEquals(expected, sum, DELTA);
                }
            }
        }
    }
}
