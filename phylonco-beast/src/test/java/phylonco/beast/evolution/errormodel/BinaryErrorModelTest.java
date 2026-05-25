package phylonco.beast.evolution.errormodel;

import beast.base.evolution.datatype.Binary;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.RealScalarParam;
import org.junit.Test;

import static junit.framework.Assert.assertEquals;

public class BinaryErrorModelTest {

    private static double DELTA = 1e-10;

    @Test
    public void testBinaryErrorModelSumsToOne() {
        Binary datatype = new Binary();

        BinaryErrorModel errorModel = new BinaryErrorModel();
        errorModel.initByName(
                "alpha", new RealScalarParam(0.1, UnitInterval.INSTANCE),
                "beta", new RealScalarParam(0.2, UnitInterval.INSTANCE),
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

    @Test
    public void testBinaryAmbiguitiesEqualOne() {
        Binary datatype = new Binary();

        BinaryErrorModel errorModel = new BinaryErrorModel();
        errorModel.initByName(
                "alpha", new RealScalarParam(0.1, UnitInterval.INSTANCE),
                "beta", new RealScalarParam(0.2, UnitInterval.INSTANCE),
                "datatype", datatype
        );
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
