package phylonco.lphy.evolution.datatype;

import jebl.evolution.sequences.Nucleotides;
import jebl.evolution.sequences.State;
import lphy.base.evolution.datatype.DataType;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * For unphased genotype data.
 * @author Alexei Drummond
 * @author Kylie Chen
 * @author Walter Xie
 */
public class UnphasedGenotype extends DataType {

    public static final String NAME = "nucleotideDiploid10";
    public static final int CANONICAL_STATE_COUNT = 10;
    public static final int STATE_COUNT = 12;

    public static final UnphasedGenotypeState[] CANONICAL_STATES;
    public static final UnphasedGenotypeState[] STATES;
    public static final UnphasedGenotypeState UNKNOWN_STATE;
    public static final UnphasedGenotypeState GAP_STATE;

    static {
        CANONICAL_STATES = new UnphasedGenotypeState[CANONICAL_STATE_COUNT];
        String codeString = "AMRWCSYGKT";

        int x = 0;
        char code;
        for(int i = 0; i < 4; i++) {
            for(int j = i; j < 4; j++) {
                String name = "" + DataType.NUCL_CHAR[i] + DataType.NUCL_CHAR[j];
                //code = (char) (x + '0');
                code = codeString.charAt(x);
                CANONICAL_STATES[x] = new UnphasedGenotypeState(name, Character.toString(code), x);
                x++;
            }
        }
        assert x == CANONICAL_STATE_COUNT;

        // no ambiguous states

        UNKNOWN_STATE = new UnphasedGenotypeState("unknown genotype", "?", 10, CANONICAL_STATES);
        GAP_STATE = new UnphasedGenotypeState("gap", "-", 11, CANONICAL_STATES);
        STATES = new UnphasedGenotypeState[STATE_COUNT];

        int i;
        for(i = 0; i < CANONICAL_STATE_COUNT; ++i) {
            STATES[i] = CANONICAL_STATES[i];
        }

        STATES[10] = UNKNOWN_STATE;
        STATES[11] = GAP_STATE;

    }

    /**
     *
     * @param parent1_index must be the first parent canonical state index
     * @param parent2_index must be the second parent canonical state index
     * @return the phasedGenotype state
     */
    public static int getUnphasedGenotypeIndex(int parent1_index, int parent2_index) {
        if (parent1_index < 0 || parent1_index >= 4 ||
                parent2_index < 0 || parent2_index >= 4) {
            throw new RuntimeException("The parents should be canonical states.");
        }
        int a = Math.min(parent1_index, parent2_index);
        int b = Math.max(parent1_index, parent2_index);

        switch (a) {
            case 0:
                return b;              // AA AC AG AT -> 0 1 2 3
            case 1:
                return 4 + (b - 1);    // CC CG CT -> 4 5 6
            case 2:
                return 7 + (b - 2);    // GG GT -> 7 8
            case 3:
                return 9;              // TT -> 9
            default:
                throw new AssertionError("The parents should be canonical states.");
        }
    }

    /**
     *
     * @param stateIndex the phased genotype state index
     * @return The two indices array for the two parents nucleotide states. The first index is the fist parent index and the
     * second index is the second parent index.
     */
    public static int[] getNucleotideIndex(int stateIndex) {

        switch (stateIndex) {
            case 0: return new int[]{0, 0}; // AA
            case 1: return new int[]{0, 1}; // AC
            case 2: return new int[]{0, 2}; // AG
            case 3: return new int[]{0, 3}; // AT
            case 4: return new int[]{1, 1}; // CC
            case 5: return new int[]{1, 2}; // CG
            case 6: return new int[]{1, 3}; // CT
            case 7: return new int[]{2, 2}; // GG
            case 8: return new int[]{2, 3}; // GT
            case 9: return new int[]{3, 3}; // TT

            case 10:
                int r = Nucleotides.getState(
                        UnphasedGenotype.INSTANCE.getState(stateIndex).getCode()
                ).getIndex();
                return new int[]{r, r};

            case 11:
                int y = Nucleotides.getState(
                        UnphasedGenotype.INSTANCE.getState(stateIndex).getCode()
                ).getIndex();
                return new int[]{y, y};

            // ...

            case 16:
                int unknown = Nucleotides.getUnknownState().getIndex();
                return new int[]{unknown, unknown};

            case 17:
                int gap = Nucleotides.getGapState().getIndex();
                return new int[]{gap, gap};

            default:
                throw new IllegalArgumentException(
                        "The unphased genotype state index is out of range.");
        }
    }

    //*** Singleton ***//

    public static UnphasedGenotype INSTANCE = new UnphasedGenotype();
    private UnphasedGenotype(){}

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
        return "GT10"; // trimmed in studio if too long
    }

    /**
     *
     * @param index index of state
     * @return canonical state at index
     */
    public static UnphasedGenotypeState getCanonicalState(int index) {
        return CANONICAL_STATES[index];
    }
}
