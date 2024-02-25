package phylonco.lphy.evolution.datatype;


import jebl.evolution.sequences.NucleotideState;
import jebl.evolution.sequences.Nucleotides;

import static jebl.evolution.sequences.Nucleotides.getState;

public class NucleotideGenotypeHelper {

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
            parent1_index = getState(originalCode).getIndex();
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
            int index = 4 * parent1_index + parent2_index;
            return index;
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
        NucleotideState[] parents_state = {getState(parent1_index),getState(parent2_index)};

        return parents_state;
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
        PhasedGenotypeState phasedState = (PhasedGenotypeState) PhasedGenotype.INSTANCE.getState(phasedIndex);

        return phasedState;
    }
}
