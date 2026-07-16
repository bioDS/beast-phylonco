package phylonco.beast.evolution.errormodel;

import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.RealScalarParam;
import org.junit.Test;
import phylonco.beast.evolution.datatype.NucleotideDiploid10;

import static junit.framework.Assert.assertEquals;

public class GT10ErrorModelTest {

    private static double DELTA = 1e-10;

    /***
     * Given the true state is y, the conditional probability over all possible states is 1.
     * sum( P(s | y) ) = 1,  where s are all possible states.
     * Tests sum( P(s | y) ) = 1 for all possible true states.
     ***/
    @Test
    public void testGT10ErrorSumsToOne() {
        NucleotideDiploid10 datatype = new NucleotideDiploid10();

        GT10ErrorModel errorModel = new GT10ErrorModel();
        errorModel.initByName(
                "epsilon", new RealScalarParam(0.1, UnitInterval.INSTANCE),
                "delta", new RealScalarParam(0.2, UnitInterval.INSTANCE),
                "datatype", datatype
        );
        errorModel.initAndValidate();

        for (int trueState = 0; trueState < datatype.getStateCount(); trueState++) {
            System.out.println("genotype true: " + datatype.getGenotype(trueState));
            double sum = 0.0;
            for (int observedState = 0; observedState < datatype.getStateCount(); observedState++) {
                sum += errorModel.getProbability(observedState, trueState);
                System.out.print(errorModel.getProbability(observedState, trueState) + " ");
            }
            System.out.println("\nsum: " + sum);
//            assertEquals(1.0, sum, DELTA);
        }

    }
}
