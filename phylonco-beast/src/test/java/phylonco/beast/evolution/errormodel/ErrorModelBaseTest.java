package phylonco.beast.evolution.errormodel;

import beast.base.evolution.datatype.Nucleotide;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.RealScalarParam;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class ErrorModelBaseTest {

    private static double DELTA = 1e-10;


    @Test
    public void testNucleotideErrorModelSumsToOne() {
        Nucleotide datatype = new Nucleotide();

        ErrorModelBase errorModel = new ErrorModelBase();
        errorModel.initByName(
                "epsilon", new RealScalarParam(0.1, UnitInterval.INSTANCE),
                "datatype", datatype);
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
        errorModel.initByName(
                "epsilon", new RealScalarParam(0.1, UnitInterval.INSTANCE),
                "datatype", datatype);
        errorModel.initAndValidate();

        for (int observedState = 0; observedState < datatype.mapCodeToStateSet.length; observedState++) {
            // observed state is ambiguous
            if (datatype.isAmbiguousCode(observedState)) {
                // all ambiguous states have probability 1.0
                int[] stateSet = datatype.getStatesForCode(observedState);
                for (int s: stateSet) {
                    double expected = 1.0;
                    double delta = 0.0;
                    double modelProb = errorModel.getProbability(observedState, s);
                    assertEquals(expected, modelProb, delta);
                }
            }
        }

    }
}

