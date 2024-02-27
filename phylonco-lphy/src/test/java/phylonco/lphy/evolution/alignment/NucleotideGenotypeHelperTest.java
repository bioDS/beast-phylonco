package phylonco.lphy.evolution.alignment;

import jebl.evolution.sequences.NucleotideState;
import jebl.evolution.sequences.State;
import org.junit.Test;
import phylonco.lphy.evolution.datatype.PhasedGenotypeState;

import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;

import static org.junit.Assert.*;
import static phylonco.lphy.evolution.datatype.PhasedGenotype.*;

public class NucleotideGenotypeHelperTest {
    /**
     * Test getNucleotideIndex function in {@link phylonco.lphy.evolution.datatype.PhasedGenotype}.
     *
     * @Paremeters phased genotype state index
     * @Returns the array of two parent indices
     */
    @Test
    public void getNucleotideIndexTest() {
        int stateIndex = 0;
        for (int i = 0; i < 24; i++) {
            int[] observed = getNucleotideIndex(stateIndex);
            if (stateIndex == 0) {
                int[] expected = {0, 0};
                assertArrayEquals("StateIndex is 0", expected, observed);
            } else if (stateIndex == 1) {
                int[] expected = {0, 1};
                assertArrayEquals("StateIndex is 1", expected, observed);
            } else if (stateIndex == 2) {
                int[] expected = {0, 2};
                assertArrayEquals("StateIndex is 2", expected, observed);
            } else if (stateIndex == 3) {
                int[] expected = {0, 3};
                assertArrayEquals("StateIndex is 3", expected, observed);
            } else if (stateIndex == 4) {
                int[] expected = {1, 0};
                assertArrayEquals("StateIndex is 4", expected, observed);
            } else if (stateIndex == 5) {
                int[] expected = {1, 1};
                assertArrayEquals("StateIndex is 5", expected, observed);
            } else if (stateIndex == 6) {
                int[] expected = {1, 2};
                assertArrayEquals("StateIndex is 6", expected, observed);
            } else if (stateIndex == 7) {
                int[] expected = {1, 3};
                assertArrayEquals("StateIndex is 7", expected, observed);
            } else if (stateIndex == 8) {
                int[] expected = {2, 0};
                assertArrayEquals("StateIndex is 8", expected, observed);
            } else if (stateIndex == 9) {
                int[] expected = {2, 1};
                assertArrayEquals("StateIndex is 9", expected, observed);
            } else if (stateIndex == 10) {
                int[] expected = {2, 2};
                assertArrayEquals("StateIndex is 10", expected, observed);
            } else if (stateIndex == 11) {
                int[] expected = {2, 3};
                assertArrayEquals("StateIndex is 11", expected, observed);
            } else if (stateIndex == 12) {
                int[] expected = {3, 0};
                assertArrayEquals("StateIndex is 12", expected, observed);
            } else if (stateIndex == 13) {
                int[] expected = {3, 1};
                assertArrayEquals("StateIndex is 13", expected, observed);
            } else if (stateIndex == 14) {
                int[] expected = {3, 2};
                assertArrayEquals("StateIndex is 14", expected, observed);
            } else if (stateIndex == 15) {
                int[] expected = {3, 3};
                assertArrayEquals("StateIndex is 15", expected, observed);
            } else if (stateIndex == 16) {
                int[] expected = {6, 6};
                assertArrayEquals("StateIndex is 16", expected, observed);
            } else if (stateIndex == 17) {
                int[] expected = {4, 4};
                assertArrayEquals("StateIndex is 17", expected, observed);
            } else if (stateIndex == 18) {
                int[] expected = {7, 7};
                assertArrayEquals("StateIndex is 17", expected, observed);
            } else if (stateIndex == 19) {
                int[] expected = {8, 8};
                assertArrayEquals("StateIndex is 19", expected, observed);
            } else if (stateIndex == 20) {
                int[] expected = {5, 5};
                assertArrayEquals("StateIndex is 20", expected, observed);
            } else if (stateIndex == 21) {
                int[] expected = {9, 9};
                assertArrayEquals("StateIndex is 21", expected, observed);
            } else if (stateIndex == 22) {
                int[] expected = {15, 15};
                assertArrayEquals("StateIndex is 22", expected, observed);
            } else if (stateIndex == 23) {
                int[] expected = {16, 16};
                assertArrayEquals("StateIndex is 23", expected, observed);
            }
            stateIndex++;
            if (stateIndex == -1) {
                String expectedErrorMessage = "The phased genotype state index should be in the range of 0 to 23!";
                try {
                    // use the function and expect to get the exception
                    getNucleotideIndex(stateIndex);

                    // the test would fail if there is no exception been thrown
                    fail("Expected an exception to be thrown");
                } catch (Exception e) {
                    assertEquals("Exception message", expectedErrorMessage, e.getMessage());
                }
            }
        }
    }

    /**
     * Test getPhasedGenotypeIndex function in {@link phylonco.lphy.evolution.datatype.PhasedGenotype}
     *
     * @Parameters first and second parent indices
     * @Returns phased genotype state index
     */
    @Test
    public void getPhasedGenotypeIndexTest() {
        int parent1_index;
        int parent2_index;
        for (parent1_index = 0; parent1_index < 4; parent1_index++) {
            for (parent2_index = 0; parent2_index < 4; parent2_index++) {
                int observed = getPhasedGenotypeIndex(parent1_index, parent2_index);

                if (parent1_index == 0) {
                    if (parent2_index == 0) {
                        int expected = 0;
                        assertEquals("0-0", expected, observed);
                    } else if (parent2_index == 1) {
                        int expected = 1;
                        assertEquals("0-1", expected, observed);
                    } else if (parent2_index == 2) {
                        int expected = 2;
                        assertEquals("0-2", expected, observed);
                    } else if (parent2_index == 3) {
                        int expected = 3;
                        assertEquals("0-3", expected, observed);
                    }
                } else if (parent1_index == 1) {
                    if (parent2_index == 0) {
                        int expected = 4;
                        assertEquals("1-0", expected, observed);
                    } else if (parent2_index == 1) {
                        int expected = 5;
                        assertEquals("1-1", expected, observed);
                    } else if (parent2_index == 2) {
                        int expected = 6;
                        assertEquals("1-2", expected, observed);
                    } else if (parent2_index == 3) {
                        int expected = 7;
                        assertEquals("1-3", expected, observed);
                    }
                } else if (parent1_index == 2) {
                    if (parent2_index == 0) {
                        int expected = 8;
                        assertEquals("2-0", expected, observed);
                    } else if (parent2_index == 1) {
                        int expected = 9;
                        assertEquals("2-1", expected, observed);
                    } else if (parent2_index == 2) {
                        int expected = 10;
                        assertEquals("2-2", expected, observed);
                    } else if (parent2_index == 3) {
                        int expected = 11;
                        assertEquals("2-3", expected, observed);
                    }
                } else if (parent1_index == 3) {
                    if (parent2_index == 0) {
                        int expected = 12;
                        assertEquals("3-0", expected, observed);
                    } else if (parent2_index == 1) {
                        int expected = 13;
                        assertEquals("3-1", expected, observed);
                    } else if (parent2_index == 2) {
                        int expected = 14;
                        assertEquals("3-2", expected, observed);
                    } else if (parent2_index == 3) {
                        int expected = 15;
                        assertEquals("3-3", expected, observed);
                    }
                }
            }
        }

        if (parent1_index == -1) {
            String expectedErrorMessage = "The parents should be canonical states.";
            try {
                // use the function and expect to get the exception
                getPhasedGenotypeIndex(parent1_index, 3);

                // the test would fail if there is no exception been thrown
                fail("Expected an exception to be thrown");
            } catch (Exception e) {
                assertEquals("Exception message", expectedErrorMessage, e.getMessage());
            }
        }
    }

    /**
     * Test the function getNucleotideState in {@link phylonco.lphy.evolution.datatype.PhasedGenotype}
     *
     * @Parameters the state code of phased genotype
     * @Return the state code array of two parents
     */
    @Test
    public void getNucleotideStateTest() throws NoSuchMethodException, InvocationTargetException, InstantiationException, IllegalAccessException {
        // create a constructor to use the final class methods
        Constructor<NucleotideState> constructor = NucleotideState.class.getDeclaredConstructor(String.class, String.class, int.class, byte.class);
        constructor.setAccessible(true);
        // create a new object for phasedGenotypeState
        PhasedGenotypeState phasedGenotypeState = new PhasedGenotypeState("name", "code", -1, new State[0]);

        // start to test
        if (phasedGenotypeState.getIndex() == 0) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Adenine", "A", 0, 0b0001);
            NucleotideState nucleotideState2 = constructor.newInstance("Adenine", "A", 0, 0b0001);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 1) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Adenine", "A", 0, 0b0001);
            NucleotideState nucleotideState2 = constructor.newInstance("Cytosine", "C", 1, 0b0010);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 2) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Adenine", "A", 0, 0b0001);
            NucleotideState nucleotideState2 = constructor.newInstance("Guanine", "G", 2, 0b0100);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 3) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Adenine", "A", 0, 0b0001);
            NucleotideState nucleotideState2 = constructor.newInstance("Thymine", "T", 3, 0b1000);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 4) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Cytosine", "C", 1, 0b0010);
            NucleotideState nucleotideState2 = constructor.newInstance("Adenine", "A", 0, 0b0001);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 5) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Cytosine", "C", 1, 0b0010);
            NucleotideState nucleotideState2 = constructor.newInstance("Cytosine", "C", 1, 0b0010);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 6) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Cytosine", "C", 1, 0b0010);
            NucleotideState nucleotideState2 = constructor.newInstance("Guanine", "G", 2, 0b0100);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 7) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Cytosine", "C", 1, 0b0010);
            NucleotideState nucleotideState2 = constructor.newInstance("Thymine", "T", 3, 0b1000);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 8) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Guanine", "G", 2, 0b0100);
            NucleotideState nucleotideState2 = constructor.newInstance("Adenine", "A", 0, 0b0001);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 9) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Guanine", "G", 2, 0b0100);
            NucleotideState nucleotideState2 = constructor.newInstance("Cytosine", "C", 1, 0b0010);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 10) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Guanine", "G", 2, 0b0100);
            NucleotideState nucleotideState2 = constructor.newInstance("Guanine", "G", 2, 0b0100);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 11) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Guanine", "G", 2, 0b0100);
            NucleotideState nucleotideState2 = constructor.newInstance("Thymine", "T", 3, 0b1000);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 12) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Thymine", "T", 3, 0b1000);
            NucleotideState nucleotideState2 = constructor.newInstance("Adenine", "A", 0, 0b0001);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 13) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Thymine", "T", 3, 0b1000);
            NucleotideState nucleotideState2 = constructor.newInstance("Cytosine", "C", 1, 0b0010);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 14) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Thymine", "T", 3, 0b1000);
            NucleotideState nucleotideState2 = constructor.newInstance("Guanine", "G", 2, 0b0100);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 15) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Thymine", "T", 3, 0b1000);
            NucleotideState nucleotideState2 = constructor.newInstance("Thymine", "T", 3, 0b1000);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 16) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Adenine", "A", 0, 0b0001);
            NucleotideState nucleotideState2 = constructor.newInstance("Cytosine", "C", 1, 0b0010);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 17) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Adenine", "A", 0, 0b0001);
            NucleotideState nucleotideState2 = constructor.newInstance("Guanine", "G", 2, 0b0100);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 18) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Adenine", "A", 0, 0b0001);
            NucleotideState nucleotideState2 = constructor.newInstance("Thymine", "T", 3, 0b1000);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 19) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Cytosine", "C", 1, 0b0010);
            NucleotideState nucleotideState2 = constructor.newInstance("Guanine", "G", 2, 0b0100);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 20) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Cytosine", "C", 1, 0b0010);
            NucleotideState nucleotideState2 = constructor.newInstance("Thymine", "T", 3, 0b1000);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 21) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Guanine", "G", 2, 0b0100);
            NucleotideState nucleotideState2 = constructor.newInstance("Thymine", "T", 3, 0b1000);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 22) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("unkown state", "?", 0, 0b0010);
            NucleotideState nucleotideState2 = constructor.newInstance("Thymine", "T", 3, 0b1000);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        } else if (phasedGenotypeState.getIndex() == 23) {
            NucleotideState[] observed = getNucleotideState(phasedGenotypeState);
            NucleotideState nucleotideState1 = constructor.newInstance("Cytosine", "C", 1, 0b0100);
            NucleotideState nucleotideState2 = constructor.newInstance("Thymine", "T", 3, 0b1000);
            NucleotideState[] expected = {nucleotideState1, nucleotideState2};
            assertArrayEquals(expected, observed);
        }
    }

    /**
     * Test the function getPhasedGenotypeState in {@link phylonco.lphy.evolution.datatype.PhasedGenotype}
     *
     * @Parameters the state codes of two nucleotide parents
     * @Return the state code of phased genotype
     */
    @Test
    public void getPhasedGenotypeStateTest() throws NoSuchMethodException, InvocationTargetException, InstantiationException, IllegalAccessException {
        // create a constructor to use the final class methods
        Constructor<NucleotideState> constructor = NucleotideState.class.getDeclaredConstructor(String.class, String.class, int.class, byte.class);
        constructor.setAccessible(true);
        // create a new object for phasedGenotypeState
        PhasedGenotypeState phasedGenotypeState = new PhasedGenotypeState("name", "code", -1, CANONICAL_STATES);

        // initialise the parent indices
        NucleotideState parent1_state = constructor.newInstance("name", "code", -1, (byte) 0b11111111);
        NucleotideState parent2_state = constructor.newInstance("name", "code", -1, (byte) 0b11111111);

        // start to test
        if (parent1_state.getIndex() == 0) {
            if (parent2_state.getIndex() == 0) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("AA", "AA", 0, CANONICAL_STATES);
                assertEquals(expected, observed);
            } else if (parent2_state.getIndex() == 1) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("AC", "AC", 1, CANONICAL_STATES);
                assertEquals(expected, observed);
            } else if (parent2_state.getIndex() == 2) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("AG", "AG", 2, CANONICAL_STATES);
                assertEquals(expected, observed);
            } else if (parent2_state.getIndex() == 3) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("AT", "AT", 3, CANONICAL_STATES);
                assertEquals(expected, observed);
            }
        } else if (parent1_state.getIndex() == 1) {
            if (parent2_state.getIndex() == 0) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("CA", "CA", 4, CANONICAL_STATES);
                assertEquals(expected, observed);
            } else if (parent2_state.getIndex() == 1) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("CC", "CC", 5, CANONICAL_STATES);
                assertEquals(expected, observed);
            } else if (parent2_state.getIndex() == 2) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("CG", "CG", 6, CANONICAL_STATES);
                assertEquals(expected, observed);
            } else if (parent2_state.getIndex() == 3) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("CT", "CT", 7, CANONICAL_STATES);
                assertEquals(expected, observed);
            }
        } else if (parent1_state.getIndex() == 2) {
            if (parent2_state.getIndex() == 0) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state); // wrong
                PhasedGenotypeState expected = new PhasedGenotypeState("GC", "GC", 8, CANONICAL_STATES);
                assertEquals(expected, observed);
            } else if (parent2_state.getIndex() == 1) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("GC", "GC", 9, CANONICAL_STATES);
                assertEquals(expected, observed);
            } else if (parent2_state.getIndex() == 2) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("GG", "GG", 10, CANONICAL_STATES);
                assertEquals(expected, observed);
            } else if (parent2_state.getIndex() == 3) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("GT", "GT", 11, CANONICAL_STATES);
                assertEquals(expected, observed);
            }
        } else if (parent1_state.getIndex() == 3) {
            if (parent2_state.getIndex() == 0) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("TA", "TA", 12, CANONICAL_STATES);
                assertEquals(expected, observed);
            } else if (parent2_state.getIndex() == 1) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("TC", "TC", 13, CANONICAL_STATES);
                assertEquals(expected, observed);
            } else if (parent2_state.getIndex() == 2) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("TG", "TG", 14, CANONICAL_STATES);
                assertEquals(expected, observed);
            } else if (parent2_state.getIndex() == 3) {
                PhasedGenotypeState observed = getPhasedGenotypeState(parent1_state, parent2_state);
                PhasedGenotypeState expected = new PhasedGenotypeState("TT", "TT", 15, CANONICAL_STATES);
                assertEquals(expected, observed);
            }
        }
    }
}