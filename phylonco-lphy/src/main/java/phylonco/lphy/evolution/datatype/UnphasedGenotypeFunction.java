package phylonco.lphy.evolution.datatype;

import jebl.evolution.sequences.SequenceType;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;

public class UnphasedGenotypeFunction extends DeterministicFunction<SequenceType> {

    public UnphasedGenotypeFunction() {}

    @GeneratorInfo(name = "unphasedGenotype",
            verbClause = "is",
            narrativeName = "unphased genotype data type",
            description = "The unphased genotype data type.")
    public Value<SequenceType> apply() {
        return new Value<>(null, UnphasedGenotype.INSTANCE, this);
    }
}
