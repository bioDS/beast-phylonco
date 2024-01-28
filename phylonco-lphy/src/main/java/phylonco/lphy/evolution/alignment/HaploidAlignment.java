package phylonco.lphy.evolution.alignment;

import jebl.evolution.sequences.SequenceType;
import lphy.base.evolution.alignment.AbstractAlignment;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.base.function.io.ReaderConst;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

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

//        // get a new object for haploid taxon number
//        int haploidTaxaNum = 2*numTaxa;

        // get the names array for new alignment
        List<String> TaxaNames = new ArrayList();
        for (int k = 0; k<numTaxa ; k++){
            String names = originalAlignment.getTaxonName(k);
            TaxaNames.add(names + "-1");
            TaxaNames.add(names + "-2");
        }

        // initialise the new alignment
        Alignment newAlignment = new SimpleAlignment((Map<String, Integer>) TaxaNames,
                originalAlignment.nchar(),SequenceType.NUCLEOTIDE);

        // map the new alignment
        for (int i = 0; i<numTaxa; i++){
            // obtain the name of each site
           List<Object> sequence = new ArrayList<>();

            for (int k = 0; k < originalAlignment.nchar(); k++){
                sequence.add(originalAlignment.getState(i,k));
            }

            for (int j = 0; j<sequence.size();j++){
                // split the diploid into haploids
                if (j%2 == 0){
                    // get the index 0 to the first haploid
                    newAlignment.setState(i*2,j, sequence.toArray.getCanonicSates(j));
                } else {
                    // get the index 1 to the second haploid
                    newAlignment.setState(i*2 + 1,j,sequence.toArray.getCanonicSates(j));
                }
            }
        }

        // return to the new alignment
        return new Value<>(null, newAlignment, this);
    }
}

