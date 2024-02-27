package phylonco.lphy.evolution.alignment;

import jebl.evolution.sequences.Nucleotides;
import jebl.evolution.sequences.SequenceType;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.AbstractAlignment;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.base.function.io.ReaderConst;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import phylonco.lphy.evolution.datatype.PhasedGenotype;

import java.util.ArrayList;
import java.util.List;

import static phylonco.lphy.evolution.datatype.PhasedGenotype.getNucleotideIndex;


public class HaploidAlignment extends DeterministicFunction<Alignment> {

    // ensure the alignment is not null
    public HaploidAlignment(@ParameterInfo(name = ReaderConst.ALIGNMENT,
    description = "the haploid alignments written from phased genotype alignment" )
                                Value<AbstractAlignment> alignmentValue){
        if (alignmentValue == null) throw new IllegalArgumentException("The alignment can't be null!");
        if (! PhasedGenotype.NAME.equals(alignmentValue.value().getSequenceType().getName()))
            throw new IllegalArgumentException("Must be phased genotype alignment!");
        setParam(ReaderConst.ALIGNMENT, alignmentValue);
    }

    // for the function
    @GeneratorInfo(name = "haploid", description = "Split the phased genotype alignment into two haploid alignments.")
    @Override
    public Value<Alignment> apply() {
        // get the values for originAlignment
        Alignment originalAlignment = ((Value<Alignment>) getParams().get(ReaderConst.ALIGNMENT)).value();

        // obtain the number of taxon in the alignment
        int numTaxa = originalAlignment.ntaxa();

        // get the names array for new alignment
        List<String> TaxaNames = new ArrayList<>();
        for (int k = 0; k<numTaxa ; k++){
            String names = originalAlignment.getTaxonName(k);
            TaxaNames.add(names + "-1");
            TaxaNames.add(names + "-2");
        }
        String[] taxaNames = TaxaNames.toArray(new String[0]);

        // initialise the new alignment
        Alignment newAlignment = new SimpleAlignment(Taxa.createTaxa(taxaNames),
                originalAlignment.nchar(),SequenceType.NUCLEOTIDE);

        // map the new alignment
        for (int i = 0; i<numTaxa; i++){
            // add each states to the seq
            for (int j = 0; j < originalAlignment.nchar(); j++){
                // get the state index of each site
                int stateIndex = originalAlignment.getState(i,j);

                // convert the phased genotype into two parents nucleotides
                int[] parentsIndex = getNucleotideIndex(stateIndex);
                int parent1_index = parentsIndex[0];
                int parent2_index = parentsIndex[1];

                // map the nucleotide states into the new alignment
                newAlignment.setState(i*2, j, parent1_index);
                newAlignment.setState(i*2+1, j, parent2_index);
            }
        }

        // return to the new alignment
        return new Value<>(null, newAlignment, this);
    }
}

