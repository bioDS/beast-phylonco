package phylonco.lphy.evolution.alignment;

import jebl.evolution.sequences.SequenceType;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.core.model.Value;
import org.junit.jupiter.api.Test;
import phylonco.lphy.evolution.datatype.PhasedGenotype;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;

public class HeterozygousTest {
    // test nucleotide
    @Test
    void HeterozygousMutationTest() {
        Alignment alignment = new SimpleAlignment(Taxa.createTaxa(1), 4, SequenceType.NUCLEOTIDE);
        Value<Alignment> alignmentValue = new Value<>("alignment", alignment);
        HeterozygousMutateAlignment instance = new HeterozygousMutateAlignment(alignmentValue, new Value<>("id", 2), null);
        String sequence = alignment.getSequence(0);

        Alignment instanceAlignment = instance.sample().value();
        // check length
        assertEquals(sequence.length(), instanceAlignment.getSequence(0).length());

        // check genotype
        int diff = 0;
        for (int i = 0; i < sequence.length(); i++) {
            String originalGenotype = String.valueOf(sequence.charAt(i));
            int newGenotype = instanceAlignment.getState(0,i);
            if (originalGenotype.equals("A") && newGenotype == 0) {
                continue;
            } else if (originalGenotype.equals("C") && newGenotype == 5) {
                continue;
            } else if (originalGenotype.equals("G") && newGenotype == 10) {
                continue;
            } else if (originalGenotype.equals("T") && newGenotype == 15) {
                continue;
            } else {
                diff ++;
            }
        }
        assertEquals(2,diff);

    }

    // test diploid
    @Test
    void testPhasedGenotypes() {
        PhasedGenotype dataType = PhasedGenotype.INSTANCE;
        Alignment alignment = new SimpleAlignment(Taxa.createTaxa(1), 4, dataType);
        Value<Alignment> alignmentValue = new Value<>("alignment", alignment);
        Value<Integer> nValue = new Value<>("id", 2);
        HeterozygousMutateAlignment instance = new HeterozygousMutateAlignment(alignmentValue, nValue, null);

        Alignment instanceAlignment = instance.sample().value();
        // check length
        assertEquals(alignment.getSequence(0).length(), instanceAlignment.getSequence(0).length());

        // check genotype
        int diff = 0;
        for (int i = 0; i < alignment.nchar(); i++) {
            int originalGenotype = alignment.getState(0, i);
            int newGenotype = instanceAlignment.getState(0, i);
            if (originalGenotype == newGenotype) {
                continue;
            } else {
                diff++;
            }
        }

        assertEquals(2, diff);
    }


    @Test
    void testSampleCanonicalState() {
        int[] ref = new int[]{0};
        int alt = HeterozygousMutateAlignment.getRandomCanonicalState(ref);

        assertNotEquals(ref[0], alt);

    }

    @Test
    void testTwoElementsForCanonicalState() {
        int[] ref = new int[]{0, 1};
        int alt = HeterozygousMutateAlignment.getRandomCanonicalState(ref);

        assertNotEquals(ref[0], alt);
        assertNotEquals(ref[1], alt);

    }

    // test with positions
    @Test
    void positionTest() {
        PhasedGenotype dataType = PhasedGenotype.INSTANCE;
        Alignment alignment = new SimpleAlignment(Taxa.createTaxa(1), 6, dataType);
        Value<Alignment> alignmentValue = new Value<>("alignment", alignment);
        Value<Integer> nValue = new Value<>("id", 4);
        Integer[] positions = new Integer[2];
        positions[0] = 0;
        positions[1] = 3;
        Value<Integer[]> positionValue = new Value<>("id", positions);
        HeterozygousMutateAlignment instance = new HeterozygousMutateAlignment(alignmentValue, nValue, null);

        Alignment instanceAlignment = instance.sample().value();
        // check length
        assertEquals(alignment.getSequence(0).length(), instanceAlignment.getSequence(0).length());

        // check genotype
        int diff = 0;
        for (int i = 0; i < alignment.nchar(); i++) {
            int originalGenotype = alignment.getState(0, i);
            int newGenotype = instanceAlignment.getState(0, i);
            if (originalGenotype == newGenotype) {
                continue;
            } else {
                diff++;
            }

            if (i == 0){
                assertNotEquals(originalGenotype, newGenotype);
            } else if (i == 3){
                assertNotEquals(originalGenotype, newGenotype);
            }
        }

        assertEquals(4, diff);
    }
}
