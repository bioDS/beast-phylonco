package phylonco.lphy.evolution.alignment;

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


public class HomozygousAlignment extends DeterministicFunction<Alignment> {
    public HomozygousAlignment(@ParameterInfo(name = ReaderConst.ALIGNMENT,
            description = "the genotype alignment (homozygous) converted from diploid alignment" )
                            Value<AbstractAlignment> alignmentValue){
        if (alignmentValue == null) throw new IllegalArgumentException("The alignment can't be null!");
        setParam(ReaderConst.ALIGNMENT, alignmentValue);
    }

    // write the function
    @GeneratorInfo(name = "homozygous", description = "convert the haploid sequence to genotype sequence.")
    @Override
    public Value<Alignment> apply() {
        // get the original seq
        Alignment originalAlignment = ((Value<Alignment>) getParams().get(ReaderConst.ALIGNMENT)).value();

        // obtain the taxa names
        List<String> taxa = new ArrayList<>();
        for (int s =0; s < originalAlignment.ntaxa();s++){
            taxa.add(originalAlignment.getTaxonName(s));
        }
        String[] taxaNames = taxa.toArray(new String[0]);

        // initialise the new alignment
        Alignment genotypeAlignment = new SimpleAlignment(Taxa.createTaxa(taxaNames),
                originalAlignment.nchar(), PhasedGenotype.INSTANCE);

        // set the alignment
        for (int i = 0; i < genotypeAlignment.ntaxa(); i++) {
            for (int j = 0; j < genotypeAlignment.nchar(); j++) {
                // get the state index of each site
                int stateIndex = originalAlignment.getState(i,j);
                // convert the nucleotide states into phased genotypes
                 int index = 4*stateIndex + stateIndex;

                 // deal with exception
//                 if (stateIndex >=4 && stateIndex <= 9 || stateIndex == 15 || stateIndex == 16 ){
//                    // ambiguous states
//                    String originalCode = Nucleotides.getState(stateIndex).getCode();
//                    index = PhasedGenotype.INSTANCE.getState(originalCode).getIndex();
//                } else if (stateIndex >9 && stateIndex <15 ){
//                    // not exist in phased genotype call them unkown state
//                    // TODO
//                    index = PhasedGenotype.INSTANCE.getUnknownState().getIndex();
//                }

                // map the new alignment states
                genotypeAlignment.setState(i,j,index);
            }
        }

        return new Value <>(null, genotypeAlignment, this);
    }
}
