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
