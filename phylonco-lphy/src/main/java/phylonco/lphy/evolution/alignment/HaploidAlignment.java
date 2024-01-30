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
import java.util.Map;

public class HaploidAlignment extends DeterministicFunction<Alignment> {

    // ensure the alignment is not null
    public HaploidAlignment(@ParameterInfo(name = ReaderConst.ALIGNMENT,
    description = "the haploid alignments written from diploid alignment" )
                                Value<AbstractAlignment> alignmentValue){
        if (alignmentValue == null) throw new IllegalArgumentException("The alignment can't be null!");
        setParam(ReaderConst.ALIGNMENT, alignmentValue);
    }

    // for the function
    @GeneratorInfo(name = "haploid", description = "Split the diploid alignment into haploid alignments.")
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

        newAlignment.taxa();

        // map the new alignment
        for (int i = 0; i<numTaxa; i++){
            // get the sequence to store the states of each seq
           List<Integer> sequence = new ArrayList<>();
            // add each states to the seq
            for (int j = 0; j < originalAlignment.nchar(); j++){

                // get the state index of each site
                int stateIndex = originalAlignment.getState(i,j);

                // convert the phased genotype states into nucleotide states
                int parent1_index = stateIndex / 4;
                int parent2_index = stateIndex % 4;

                // deal with exceptions
                if (stateIndex > 15 && stateIndex <= 21) {
                    // ambiguous (unphased) state
                    // get the code for phased state
                    String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
                    // get the nucleotide state
                    parent1_index = Nucleotides.getState(originalCode).getIndex();
                    parent2_index = parent1_index;
                } else if (stateIndex > 21) {
                    // unkown genotype and gap
                    parent1_index = Nucleotides.getGapState().getIndex();
                    parent2_index = parent1_index;
                }

                // map the nucleotide states into the new alignment
                newAlignment.setState(i*2, j, parent1_index);
                newAlignment.setState(i*2+1, j, parent2_index);
            }
        }

        // return to the new alignment
        return new Value<>(null, newAlignment, this);
    }
}

