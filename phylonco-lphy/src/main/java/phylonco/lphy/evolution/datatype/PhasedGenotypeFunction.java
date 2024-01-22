package phylonco.lphy.evolution.datatype;

import jebl.evolution.sequences.SequenceType;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;

public class PhasedGenotypeFunction extends DeterministicFunction<SequenceType> {

    public PhasedGenotypeFunction() {}

    @GeneratorInfo(name = "phasedGenotype",
            verbClause = "is",
            narrativeName = "phased genotype data type",
            description = "The phased genotype data type.")
    public Value<SequenceType> apply() {
        return new Value<>(null, PhasedGenotype.INSTANCE, this);
    }
}
