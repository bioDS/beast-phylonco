package phylonco.lphy.io;

import lphy.base.evolution.alignment.Alignment;
import lphy.base.function.io.ReaderConst;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import phylonco.lphy.evolution.alignment.HaploidAlignment;

public class WriteHaploidFasta extends DeterministicFunction<HaploidAlignment> {
    public WriteHaploidFasta(@ParameterInfo(name = ReaderConst.ALIGNMENT,
            description = "the GT16 diploid alignment that is written into the haploid fasta file.")
                             Value<Alignment> alignmentValue ) {

        if (alignmentValue == null) throw new IllegalArgumentException("The alignment can't be null!");
        setParam(ReaderConst.ALIGNMENT, alignmentValue);
    }

    @GeneratorInfo(name="HaploidType", verbClause = "is read from", narrativeName = "fasta file",
            category = GeneratorCategory.TAXA_ALIGNMENT, examples = {"covidDPG.lphy"},
            description = "A function that split diploid alignment into two haploid alignments and generate a fasta file.")
    @Override
    public Value<HaploidAlignment> apply() {
        Alignment alignment = ((Value<Alignment>) getParams().get(ReaderConst.ALIGNMENT)).value();

        // this only creates taxa and nchar
        HaploidAlignment faData = new HaploidAlignment(alignment.nchar(), alignment);

        // fill in states
        for (int i=0; i < alignment.ntaxa(); i++) {
            for (int j = 0; j < alignment.nchar(); j++) {
                faData.setState(i, j, alignment.getState(i, j));
            }
        }

        return new Value<>(null, faData, this);
    }
}

