package phylonco.lphy.evolution.alignment;

import jebl.evolution.sequences.Nucleotides;
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

        // initialise the new alignment
        Alignment genotypeAlignment = new SimpleAlignment(Taxa.createTaxa(originalAlignment.ntaxa()),
                originalAlignment.nchar(), PhasedGenotype.INSTANCE);

        // set the alignment
        for (int i = 0; i < genotypeAlignment.ntaxa(); i++) {
            for (int j = 0; j < genotypeAlignment.nchar(); j++) {
                // get the state index of each site
                int stateIndex = originalAlignment.getState(i,j);

                // convert the nucleotide states into phased genotypes
                int index = stateIndex^2 + 3;

                // deal with exceptions
                if (stateIndex >=4 && stateIndex <= 9 && stateIndex == 15 && stateIndex == 16 ){
                    // ambiguous states
                    String originalCode = Nucleotides.getState(stateIndex).getCode();
                    index = PhasedGenotype.INSTANCE.getState(originalCode).getIndex();
                } else {
                    // not exist in phased genotype ? how to deal with this
                    index = PhasedGenotype.INSTANCE.getGapState().getIndex();
                }

                // map the new alignment states
                genotypeAlignment.setState(i,j,index);
            }
        }

        return new Value <>(null, genotypeAlignment, this);
    }

    private int homozygote (int state) {
        switch (state) {
            case 0:
                return 0; // A --> AA
            case 1:
                return 4; // C --> CC
            case 2:
                return 7; // G --> GG
            case 3:
                return 9; // T --> TT
        }
        throw new RuntimeException("Unexpected state: " + state);
    }
}
