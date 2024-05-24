package phylonco.lphy.evolution.alignment;

import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.core.model.Value;
import org.junit.Test;
import phylonco.lphy.evolution.datatype.PhasedGenotype;

import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.assertEquals;

public class UnphaseGenotypeAlignmentTest {

    @Test
    public void testUnphaseAlignment() {

        Map<String, Integer> taxa = new HashMap<>();
        taxa.put("taxon1", 0);
        taxa.put("taxon2", 1);
        int nchar = 2;
        PhasedGenotype dataType = PhasedGenotype.INSTANCE;

        SimpleAlignment mockAlignment = new SimpleAlignment(taxa, nchar, dataType);

        mockAlignment.setState(0, 0, 1);  // AC
        mockAlignment.setState(0, 1, 2);  // AG
        mockAlignment.setState(1, 0, 3);  // AT
        mockAlignment.setState(1, 1, 4);  // CA

        Value<Alignment> mockValue = new Value<>(null, mockAlignment, null);

        UnphaseGenotypeAlignment unphaseGenotypeAlignment = new UnphaseGenotypeAlignment(mockValue);

        Value<Alignment> unphasedValue = unphaseGenotypeAlignment.apply();

        Alignment unphasedAlignment = unphasedValue.value();

        // Verify the unphased states
        assertEquals(16, unphasedAlignment.getState(0, 0)); // AC -> AC/CA
        assertEquals(17, unphasedAlignment.getState(0, 1)); // AG -> AG/GA
        assertEquals(18, unphasedAlignment.getState(1, 0)); // AT -> AT/TA
        assertEquals(16, unphasedAlignment.getState(1, 1)); // CA -> AC/CA
    }

}
