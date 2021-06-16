package test.beast.evolution.errormodel;

import beast.evolution.datatype.DataType;
import beast.evolution.datatype.NucleotideDiploid16;
import beast.evolution.errormodel.GT16ErrorModel;
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
                sum += errorModel.getProbability(observedState, trueState);
            }
            assertEquals(1.0, sum, DELTA);
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
                double modelProb = errorModel.getProbability(observedState, trueState);
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
                            double modelProb = errorModel.getProbability(otherStates, trueState);
                            sum += modelProb;
                        }
                    }
                    double expected = 1.0;
                    double modelProb = errorModel.getProbability(observedState, trueState);
                    sum += modelProb;
                    assertEquals(expected, sum, DELTA);
                }
            }
        }
    }
}
