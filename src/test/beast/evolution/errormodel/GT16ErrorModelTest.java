package test.beast.evolution.errormodel;

import beast.evolution.datatype.NucleotideDiploid16;
import beast.evolution.errormodel.GT16ErrorModel;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

public class GT16ErrorModelTest {

    private static double DELTA = 1e-10;

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
}
