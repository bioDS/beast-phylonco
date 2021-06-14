package test.beast.evolution.errormodel;

import beast.evolution.datatype.Nucleotide;
import beast.evolution.errormodel.ErrorModelBase;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

public class ErrorModelBaseTest {

    private static double DELTA = 1e-10;

    @Test
    public void testNucleotideErrorModelSumsToOne() {
        Nucleotide datatype = new Nucleotide();

        ErrorModelBase errorModel = new ErrorModelBase();
        errorModel.initByName("epsilon", "0.1", "datatype", datatype);
        errorModel.initAndValidate();

        for (int trueState = 0; trueState < datatype.getStateCount(); trueState++) {
            double sum = 0.0;
            for (int observedState = 0; observedState < datatype.getStateCount(); observedState++) {
                sum += errorModel.getProbability(observedState, trueState);
            }
            assertEquals(1.0, sum, DELTA);
        }

    }

    @Test
    public void testNucleotideAmbiguitiesEqualOne() {
        Nucleotide datatype = new Nucleotide();

        ErrorModelBase errorModel = new ErrorModelBase();
        errorModel.initByName("epsilon", "0.1", "datatype", datatype);
        errorModel.initAndValidate();

        for (int observedState = 0; observedState < datatype.mapCodeToStateSet.length; observedState++) {
            // observed state is ambiguous
            if (datatype.isAmbiguousCode(observedState)) {
                // all true states have probability 1.0
                for (int trueState = 0; trueState < datatype.getStateCount(); trueState++) {
                    double expected = 1.0;
                    double delta = 0.0;
                    double modelProb = errorModel.getProbability(observedState, trueState);
                    assertEquals(expected, modelProb, delta);
                }
            }
        }

    }
}

