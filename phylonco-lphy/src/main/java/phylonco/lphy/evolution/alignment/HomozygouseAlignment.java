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


public class HomozygouseAlignment extends DeterministicFunction<Alignment> {
    public HomozygouseAlignment(@ParameterInfo(name = ReaderConst.ALIGNMENT,
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
                genotypeAlignment.setState(i, j, homozygote(originalAlignment.getState(i, j)));
            }
        }

        return new Value <>(null, genotypeAlignment, this);
    }
    private int homozygote (int state) {
        switch (state) {
            case 0:
                return 0; // A --> AA
            case 1:
                return ; // C --> CC
            case 2:
                return ; // T --> TT
            case 3:
                return ; // G --> GG
        }
        throw new RuntimeException("Unexpected state: " + state);
    }
}
