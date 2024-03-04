
package phylonco.lphy.evolution.datatype;

import jebl.evolution.sequences.NucleotideState;
import jebl.evolution.sequences.Nucleotides;
import jebl.evolution.sequences.State;
import lphy.base.evolution.datatype.DataType;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * For phased genotype data.
 * @author Alexei Drummond
 * @author Kylie Chen
 * @author Walter Xie
 */
public final class PhasedGenotype extends DataType {

    public static final String NAME = "nucleotideDiploid16";
    public static final int CANONICAL_STATE_COUNT = 16;
    public static final int STATE_COUNT = 24;

    public static final PhasedGenotypeState[] CANONICAL_STATES;
    public static final PhasedGenotypeState[] STATES;
    public static final PhasedGenotypeState UNKNOWN_STATE;
    public static final PhasedGenotypeState GAP_STATE;
    public static final PhasedGenotypeState AC_OR_CA;
    public static final PhasedGenotypeState AG_OR_GA;
    public static final PhasedGenotypeState AT_OR_TA;
    public static final PhasedGenotypeState CG_OR_GC;
    public static final PhasedGenotypeState CT_OR_TC;
    public static final PhasedGenotypeState GT_OR_TG;

    static {
        CANONICAL_STATES = new PhasedGenotypeState[CANONICAL_STATE_COUNT];

        int x = 0;
        char code;
        // change from ++i to i++ for easier reading
        for(int i = 0; i < 4; i++) {
            for(int j = 0; j < 4; j++) {
                String name = "" + DataType.NUCL_CHAR[i] + DataType.NUCL_CHAR[j];
                code = x < 10 ? (char) (x + '0') : (char) (x - 10 + 'a');
                CANONICAL_STATES[x] = new PhasedGenotypeState(name, Character.toString(code), x);
                x = x + 1;
            }
        }
        assert x == CANONICAL_STATE_COUNT;

        AC_OR_CA = new PhasedGenotypeState("ac", "M", 16, CANONICAL_STATES);
        AG_OR_GA = new PhasedGenotypeState("ag", "R", 17, CANONICAL_STATES);
        AT_OR_TA = new PhasedGenotypeState("at", "W", 18, CANONICAL_STATES);
        CG_OR_GC = new PhasedGenotypeState("cg", "S", 19, CANONICAL_STATES);
        CT_OR_TC = new PhasedGenotypeState("ct", "Y", 20, CANONICAL_STATES);
        GT_OR_TG = new PhasedGenotypeState("gt", "K", 21, CANONICAL_STATES);
        UNKNOWN_STATE = new PhasedGenotypeState("unknown genotype", "?", 22, CANONICAL_STATES);
        GAP_STATE = new PhasedGenotypeState("gap", "-", 23, CANONICAL_STATES);
        STATES = new PhasedGenotypeState[STATE_COUNT];

        int i;
        for(i = 0; i < CANONICAL_STATE_COUNT; ++i) {
            STATES[i] = CANONICAL_STATES[i];
        }

        STATES[16] = AC_OR_CA;
        STATES[17] = AG_OR_GA;
        STATES[18] = AT_OR_TA;
        STATES[19] = CG_OR_GC;
        STATES[20] = CT_OR_TC;
        STATES[21] = GT_OR_TG;
        STATES[22] = UNKNOWN_STATE;
        STATES[23] = GAP_STATE;
    }

    //*** Singleton ***//

    public static PhasedGenotype INSTANCE = new PhasedGenotype();
    private PhasedGenotype(){}

    //*** implementations ***//

    @Override
    public int getStateCount() {
        return STATE_COUNT;
    }

    @Override
    public List<State> getStates() {
        return Collections.unmodifiableList(Arrays.asList((State[])STATES));
    }

    @Override
    public int getCanonicalStateCount() {
        return CANONICAL_STATE_COUNT;
    }

    @Override
    public List<? extends State> getCanonicalStates() {
        return Collections.unmodifiableList(Arrays.asList((State[])CANONICAL_STATES));
    }

    @Override
    public int getCodeLength() {
        return 2;
    }

    @Override
    public State getState(int index) {
        return STATES[index];
    }

    @Override
    public State getUnknownState() {
        return UNKNOWN_STATE;
    }

    @Override
    public State getGapState() {
        return GAP_STATE;
    }

    /**
     * @return true if the given state represents unknown genotype state.
     */
    @Override
    public boolean isUnknown(State state) {
        return state == UNKNOWN_STATE;
    }

    @Override
    public boolean isGap(State state) {
        return state == GAP_STATE;
    }

    @Override
    public String getName() {
        return NAME;
    }

    @Override
    public String getNexusDataType() {
        return NAME;
    }

    @Override
    public String toString() {
        return "GT16"; // trimmed in studio if too long
    }

    /**
     *
     * @param stateIndex the phased genotype state index
     * @return The two indices array for the two parents nucleotide states. The first index is the fist parent index and the
     * second index is the second parent index.
     */
    public static int[] getNucleotideIndex(int stateIndex){
        // initialise the parents indices
        int parent1_index = -1;
        int parent2_index = -1;
        // deal with the state indices
        if (stateIndex>=0 && stateIndex < 16){
            parent1_index = stateIndex / 4;
            parent2_index = stateIndex % 4;
        } else if (stateIndex >= 16 && stateIndex <= 21) {
            // get the code for phased state (e.g. "R" for AG or GA)
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state for the ambiguous code
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex == 22) {
            // unkown genotype
            parent1_index = Nucleotides.getUnknownState().getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex == 23){
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        } else {
            throw new IllegalArgumentException("The phased genotype state index should be in the range of 0 to 23!");
        }
        return new int[]{parent1_index,parent2_index};
    }

    /**
     *
     * @param parent1_index must be the first parent canonical state index
     * @param parent2_index must be the second parent canonical state index
     * @return the phasedGenotype state
     */
    public static int getPhasedGenotypeIndex(int parent1_index,int parent2_index){
        if (parent1_index <4 && parent1_index>=0 && parent2_index>=0 && parent2_index<4){
            return 4 * parent1_index + parent2_index;
        }  throw new RuntimeException( "The parents should be canonical states." );
    }


    /**
     *
     * @param state the state code of phased genotype
     * @return the state codes array of two nucleotide parents
     */
    public static NucleotideState[] getNucleotideState(PhasedGenotypeState state){
        // obtain the parents index array
        int[] parentsIndex = getNucleotideIndex(state.getIndex());

        // extract each parent index
        int parent1_index = parentsIndex[0];
        int parent2_index = parentsIndex[1];

        // convert the index into state and generate the array we want
        return new NucleotideState[]{Nucleotides.getState(parent1_index),Nucleotides.getState(parent2_index)};
    }

    /**
     *
     * @param parent1 must be the first parent nucleotide canonical state code
     * @param parent2 must be the second parent nucleotide canonical state code
     * @return phased genotype state code
     */
    public static PhasedGenotypeState getPhasedGenotypeState(NucleotideState parent1, NucleotideState parent2){
        // obtain the phasedGenotype index
        int phasedIndex = getPhasedGenotypeIndex(parent1.getIndex() , parent2.getIndex());

        // convert the phased index into state
        return (PhasedGenotypeState) PhasedGenotype.INSTANCE.getState(phasedIndex);
    }
}
