package phylonco.lphy.evolution.alignment;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.core.model.Value;
import org.junit.Test;
import phylonco.lphy.evolution.datatype.UnphasedGenotype;

import static org.junit.Assert.assertEquals;

public class GT10ErrorModelTest {

    private static double DELTA = 1e-15;

    private Value<Alignment> getDummySequence() {
        int numTaxa = 10;
        int numChar = 200;
        Taxa taxa = Taxa.createTaxa(numTaxa);
        UnphasedGenotype seqType = UnphasedGenotype.INSTANCE;
        Value<Alignment> alignment = new Value<>("alignment", new SimpleAlignment(taxa, numChar, seqType));
        return alignment;
    }


    @Test
    public void errorMatrixRowSumsToOne() {
        double delta = 0.1;
        double epsilon = 0.2;

        Value<Double> ep = new Value<>("epsilon", epsilon);
        Value<Double> dt = new Value<>("delta", delta);
        GT10ErrorModel errorModel = new GT10ErrorModel(ep, dt, getDummySequence());

        double expected = 1.0;
        double[][] observedMatrix = errorModel.errorMatrix(epsilon, delta);
        for (int row = 0; row < observedMatrix.length; row++) {
            double sum = 0.0;
            for (int col = 0; col < observedMatrix[row].length; col++) {
                sum += observedMatrix[row][col];
            }
            System.out.println("row: " + row + ", sum = " + sum);
            assertEquals(expected, sum, DELTA);
        }
    }
}
