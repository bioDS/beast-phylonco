package phylonco.lphy.evolution.datatype;

import jebl.evolution.sequences.SequenceType;
import jebl.evolution.sequences.State;

public class PhasedGenotypeState extends State {

    public PhasedGenotypeState(String name, String stateCode, int index) {
        super(name, stateCode, index);
    }

    public PhasedGenotypeState(String name, String stateCode, int index, State[] ambiguities) {
        super(name, stateCode, index, ambiguities);
    }

    @Override
    public boolean isGap() {
        return this == PhasedGenotype.GAP_STATE;
    }

    @Override
    public SequenceType getType() {
        return PhasedGenotype.INSTANCE;
    }
}
