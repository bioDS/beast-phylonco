package phylonco.lphy.evolution.alignment;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.core.model.Value;
import org.junit.Test;
import phylonco.lphy.evolution.datatype.PhasedGenotype;
import phylonco.lphy.evolution.readcountmodel.MutableAlignmentModel;

import static org.junit.Assert.assertArrayEquals;

public class MutableAlignmentModelTest {
    private static double DELTA = 1e-15;


    @Test
    public void testModel() {
        int numTaxa = 10;
        int numChar = 200;
        Taxa taxa = Taxa.createTaxa(numTaxa);
        PhasedGenotype seqType = PhasedGenotype.INSTANCE;
        Value<Alignment> parent = new Value<>("alignment", new SimpleAlignment(taxa, numChar, seqType));
        MutableAlignmentModel model = new MutableAlignmentModel(parent);
        Value<Alignment> child = new Value<>("alignment", model.sample().value());
        double[] expected = new double[numTaxa * numChar];
        double[] observed = new double[numTaxa * numChar];
        for (int i = 0; i < numTaxa; i++) {
            for (int j = 0; j < numChar; j++) {
                expected[i * numChar + j] = 0;
                observed[i * numChar + j] = parent.value().getState(i,j) - child.value().getState(i,j);
            }
        }
        assertArrayEquals(expected, observed, DELTA);
    }
}
