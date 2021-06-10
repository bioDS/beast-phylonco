package test.beast.evolution.errormodel;

import beast.evolution.datatype.Binary;
import beast.evolution.errormodel.BinaryErrorModel;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

public class BinaryErrorModelTest {

    private static double DELTA = 1e-10;

    @Test
    public void testBinaryErrorModelSumsToOne() {
        Binary datatype = new Binary();

        BinaryErrorModel errorModel = new BinaryErrorModel();
        errorModel.initByName(
                "alpha", "0.1",
                "beta", "0.2",
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
