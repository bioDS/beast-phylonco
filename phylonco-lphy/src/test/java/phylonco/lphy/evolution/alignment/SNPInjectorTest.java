package phylonco.lphy.evolution.alignment;

import jebl.evolution.sequences.SequenceType;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.base.evolution.datatype.Variant;
import lphy.core.model.Value;
import org.junit.jupiter.api.Test;
import phylonco.lphy.evolution.datatype.PhasedGenotype;

import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static phylonco.lphy.evolution.datatype.PhasedGenotype.getPhasedGenotypeIndex;

public class SNPInjectorTest {
    @Test
    void test() {
        Alignment a = new SimpleAlignment(Taxa.createTaxa(1), 5, SequenceType.NUCLEOTIDE);

        Variant var1 = new Variant(a.taxaNames()[0], 2, 0, 1, "0|1");
        Variant var2 = new Variant(a.taxaNames()[0], 4, 0, 1, "1|1");

        Variant[] snps = new Variant[]{var1, var2};

        Value<Variant[]> variantsValue = new Value<>("id", snps);
        Value<Alignment> alignmentValue = new Value<>("id", a);

        SNPInjector injector = new SNPInjector(alignmentValue, variantsValue);
        Alignment alignment = injector.apply().value();

        Alignment homoAlignment = new SimpleAlignment(Taxa.createTaxa(alignment.getTaxaNames()),
                alignment.nchar(), PhasedGenotype.INSTANCE);

        for (int i = 0; i<a.nchar(); i++){
            int ref = a.getState(0,i);
            int newIndex = getPhasedGenotypeIndex(ref,ref);
            homoAlignment.setState(0,i,newIndex);
        }
        List<Variant> newVariants = new ArrayList<>();

        for (int i = 0; i<alignment.nchar();  i++){
            int oldState = homoAlignment.getState(0,i);
            int newState = alignment.getState(0,i);
            if (newState != oldState){
                int[] indices = PhasedGenotype.getNucleotideIndex(newState);
                int ref = oldState;
                int alt = -1;
                String genotype;

                if (ref == indices[0]){
                    genotype = "0|1";
                    alt = indices[1];
                } else if (ref == indices[1]) {
                    genotype = "0|1";
                    alt = indices[0];
                } else {
                    genotype = "1|1";
                }
                Variant var = new Variant(homoAlignment.taxaNames()[0], i, ref, alt , genotype);
                newVariants.add(var);
            }

        }

        assertEquals(snps.length, newVariants.size());
        assertEquals(1, alignment.getState(0,2));
        assertEquals(5, alignment.getState(0,4));

    }


}
