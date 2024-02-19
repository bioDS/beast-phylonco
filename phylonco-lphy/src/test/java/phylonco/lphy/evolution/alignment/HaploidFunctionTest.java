package phylonco.lphy.evolution.alignment;

import jebl.evolution.sequences.Nucleotides;
import lphy.base.evolution.alignment.Alignment;
import org.junit.Test;
import phylonco.lphy.evolution.datatype.PhasedGenotype;

import static org.junit.Assert.assertEquals;


public class HaploidFunctionTest {
    @Test
    public void getExpectedState0() {
        int stateIndex = 0;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For AA parent 1",0, parent1_index,0);
        assertEquals("For AA parent 2",0, parent2_index,0);
    }

    @Test
    public void getExpectedState1() {
        int stateIndex = 1;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For AC parent 1", 0,parent1_index,0);
        assertEquals("For AC parent 2",1, parent2_index,0);
    }
    @Test
    public void getExpectedState2() {
        int stateIndex = 2;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For AG parent 1", 0,parent1_index,0);
        assertEquals("For AG parent 2",2, parent2_index,0);
    }

    @Test
    public void getExpectedState3() {
        int stateIndex = 3;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For AT parent 1", 0,parent1_index,0);
        assertEquals("For AT parent 2",3, parent2_index,0);
    }

    @Test
    public void getExpectedState4() {
        int stateIndex = 4;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For CA parent 1", 1,parent1_index,0);
        assertEquals("For CA parent 2",0, parent2_index,0);
    }

    @Test
    public void getExpectedState5() {
        int stateIndex = 5;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For CC parent 1", 1,parent1_index,0);
        assertEquals("For CC parent 2",1, parent2_index,0);
    }
    @Test
    public void getExpectedState6() {
        int stateIndex = 6;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For CG parent 1", 1,parent1_index,0);
        assertEquals("For CG parent 2",2, parent2_index,0);
    }
    @Test
    public void getExpectedState7() {
        int stateIndex = 7;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For CT parent 1", 1,parent1_index,0);
        assertEquals("For CT parent 2",3, parent2_index,0);
    }

    @Test
    public void getExpectedState8() {
        int stateIndex = 8;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For GA parent 1", 2,parent1_index,0);
        assertEquals("For GA parent 2",0, parent2_index,0);
    }

    @Test
    public void getExpectedState9() {
        int stateIndex = 9;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For GC parent 1", 2,parent1_index,0);
        assertEquals("For GC parent 2",1, parent2_index,0);
    }

    @Test
    public void getExpectedState10() {
        int stateIndex = 10;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For GG parent 1", 2,parent1_index,0);
        assertEquals("For GG parent 2",2, parent2_index,0);
    }
    @Test
    public void getExpectedState11() {
        int stateIndex = 11;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For GT parent 1", 2,parent1_index,0);
        assertEquals("For GT parent 2",3, parent2_index,0);
    }
    @Test
    public void getExpectedState12() {
        int stateIndex = 12;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For TA parent 1", 3,parent1_index,0);
        assertEquals("For TA parent 2",0, parent2_index,0);
    }
    @Test
    public void getExpectedState13() {
        int stateIndex = 13;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For TC parent 1", 3,parent1_index,0);
        assertEquals("For TC parent 2",1, parent2_index,0);
    }
    @Test
    public void getExpectedState14() {
        int stateIndex = 14;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For TG parent 1", 3,parent1_index,0);
        assertEquals("For TG parent 2",2, parent2_index,0);
    }
    @Test
    public void getExpectedState15() {
        int stateIndex = 15;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For TT parent 1", 3,parent1_index,0);
        assertEquals("For TT parent 2",3, parent2_index,0);
    }
    @Test
    public void getExpectedState16() {
        int stateIndex = 16;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For M parent 1", 6,parent1_index,0);
        assertEquals("For M parent 2",6, parent2_index,0);
    }
    @Test
    public void getExpectedState17() {
        int stateIndex = 17;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For R parent 1", 4,parent1_index,0);
        assertEquals("For R parent 2",4, parent2_index,0);
    }
    @Test
    public void getExpectedState18() {
        int stateIndex = 18;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For W parent 1", 7,parent1_index,0);
        assertEquals("For W parent 2",7, parent2_index,0);
    }
    @Test
    public void getExpectedState19() {
        int stateIndex = 19;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For S parent 1", 8,parent1_index,0);
        assertEquals("For S parent 2",8, parent2_index,0);
    }
    @Test
    public void getExpectedState20() {
        int stateIndex = 20;
        // convert the phased genotype states into nucleotide states
        int parent1_index = stateIndex / 4;
        int parent2_index = stateIndex % 4;
        // deal with exceptions
        if (stateIndex > 15 && stateIndex <= 21) {
            // ambiguous (unphased) state
            // get the code for phased state
            String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
            // get the nucleotide state
            parent1_index = Nucleotides.getState(originalCode).getIndex();
            parent2_index = parent1_index;
        } else if (stateIndex > 21) {
            // unkown genotype and gap
            parent1_index = Nucleotides.getGapState().getIndex();
            parent2_index = parent1_index;
        }
        assertEquals("For Y parent 1", 5,parent1_index,0);
        assertEquals("For Y parent 2",5, parent2_index,0);
    }
    @Test
    public void getExpectedState22and23() {
        for (int stateIndex=22;stateIndex<24;stateIndex++) {
            // convert the phased genotype states into nucleotide states
            int parent1_index = stateIndex / 4;
            int parent2_index = stateIndex % 4;
            // deal with exceptions
            if (stateIndex > 15 && stateIndex <= 21) {
                // ambiguous (unphased) state
                // get the code for phased state
                String originalCode = PhasedGenotype.INSTANCE.getState(stateIndex).getCode();
                // get the nucleotide state
                parent1_index = Nucleotides.getState(originalCode).getIndex();
                parent2_index = parent1_index;
            } else if (stateIndex > 21) {
                // unkown genotype and gap
                parent1_index = Nucleotides.getGapState().getIndex();
                parent2_index = parent1_index;
            }
            assertEquals("For unkown parent 1", 16, parent1_index, 0);
            assertEquals("For unkown parent 2", 16, parent2_index, 0);
        }
    }
}

