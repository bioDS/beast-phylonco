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

        for (int i = 0; i < datatype.getStateCount(); i++) {
            double sum = 0.0;
            for (int j = 0; j < datatype.getStateCount(); j++) {
                // getProbability(observed state, true state)
                sum += errorModel.getProbability(j, i);
            }
            assertEquals(1.0, sum, DELTA);
        }

    }
}

