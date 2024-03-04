package phylonco.lphy.evolution.alignment;

import jebl.evolution.sequences.NucleotideState;
import jebl.evolution.sequences.State;
import org.junit.Test;
import phylonco.lphy.evolution.datatype.PhasedGenotype;
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

        // start to test
        PhasedGenotypeState AA = new PhasedGenotypeState("AA", "0", 0);
        NucleotideState[] observedAA = getNucleotideState(AA);
        NucleotideState nucleotideStateAA1 = constructor.newInstance("Adenine", "A", 0, (byte)0b0001);
        NucleotideState nucleotideStateAA2 = constructor.newInstance("Adenine", "A", 0, (byte)0b0001);
        NucleotideState[] expectedAA = {nucleotideStateAA1, nucleotideStateAA2};
        assertEquals(expectedAA[0].getIndex(),observedAA[0].getIndex());
        assertEquals(expectedAA[1].getIndex(),observedAA[1].getIndex());


        PhasedGenotypeState AC = new PhasedGenotypeState("AC", "1", 1);
        NucleotideState[] observedAC = getNucleotideState(AC);
        NucleotideState nucleotideStateAC1 = constructor.newInstance("Adenine", "A", 0, (byte)0b0001);
        NucleotideState nucleotideStateAC2 = constructor.newInstance("Cytosine", "C", 1, (byte)0b0010);
        NucleotideState[] expectedAC = {nucleotideStateAC1, nucleotideStateAC2};
        assertEquals(expectedAC[0].getIndex(),observedAC[0].getIndex());
        assertEquals(expectedAC[1].getIndex(),observedAC[1].getIndex());

        PhasedGenotypeState AG = new PhasedGenotypeState("AG", "2", 2);
        NucleotideState[] observedAG = getNucleotideState(AG);
        NucleotideState nucleotideStateAG1 = constructor.newInstance("Adenine", "A", 0, (byte)0b0001);
        NucleotideState nucleotideStateAG2 = constructor.newInstance("Guanine", "G", 2, (byte)0b0100);
        NucleotideState[] expectedAG = {nucleotideStateAG1, nucleotideStateAG2};
        assertEquals(expectedAG[0].getIndex(),observedAG[0].getIndex());
        assertEquals(expectedAG[1].getIndex(),observedAG[1].getIndex());

        PhasedGenotypeState AT = new PhasedGenotypeState("AA", "3", 3);
        NucleotideState[] observedAT = getNucleotideState(AT);
        NucleotideState nucleotideStateAT1 = constructor.newInstance("Adenine", "A", 0, (byte)0b0001);
        NucleotideState nucleotideStateAT2 = constructor.newInstance("Thymine", "T", 3, (byte)0b1000);
        NucleotideState[] expectedAT = {nucleotideStateAT1, nucleotideStateAT2};
        assertEquals(expectedAT[0].getIndex(),observedAT[0].getIndex());
        assertEquals(expectedAT[1].getIndex(),observedAT[1].getIndex());

        PhasedGenotypeState CA = new PhasedGenotypeState("CA", "4", 4);
        NucleotideState[] observedCA = getNucleotideState(CA);
        NucleotideState nucleotideStateCA1 = constructor.newInstance("Cytosine", "C", 1, (byte)0b0010);
        NucleotideState nucleotideStateCA2 = constructor.newInstance("Adenine", "A", 0, (byte)0b0001);
        NucleotideState[] expectedCA = {nucleotideStateCA1, nucleotideStateCA2};
        assertEquals(expectedCA[0].getIndex(),observedCA[0].getIndex());
        assertEquals(expectedCA[1].getIndex(),observedCA[1].getIndex());

        PhasedGenotypeState CC = new PhasedGenotypeState("CC", "5", 5);
        NucleotideState[] observedCC = getNucleotideState(CC);
        NucleotideState nucleotideStateCC1 = constructor.newInstance("Cytosine", "C", 1, (byte)0b0010);
        NucleotideState nucleotideStateCC2 = constructor.newInstance("Cytosine", "C", 1, (byte)0b0010);
        NucleotideState[] expectedCC = {nucleotideStateCC1, nucleotideStateCC2};
        assertEquals(expectedCC[0].getIndex(),observedCC[0].getIndex());
        assertEquals(expectedCC[1].getIndex(),observedCC[1].getIndex());

        PhasedGenotypeState CG = new PhasedGenotypeState("CG", "6", 6);
        NucleotideState[] observedCG = getNucleotideState(CG);
        NucleotideState nucleotideStateCG1 = constructor.newInstance("Cytosine", "C", 1, (byte)0b0010);
        NucleotideState nucleotideStateCG2 = constructor.newInstance("Guanine", "G", 2, (byte)0b0100);
        NucleotideState[] expectedCG = {nucleotideStateCG1, nucleotideStateCG2};
        assertEquals(expectedCG[0].getIndex(),observedCG[0].getIndex());
        assertEquals(expectedCG[1].getIndex(),observedCG[1].getIndex());


        PhasedGenotypeState CT = new PhasedGenotypeState("CT", "7", 7);
        NucleotideState[] observedCT = getNucleotideState(CT);
        NucleotideState nucleotideStateCT1 = constructor.newInstance("Cytosine", "C", 1, (byte)0b0010);
        NucleotideState nucleotideStateCT2 = constructor.newInstance("Thymine", "T", 3, (byte)0b1000);
        NucleotideState[] expectedCT = {nucleotideStateCT1, nucleotideStateCT2};
        assertEquals(expectedCT[0].getIndex(),observedCT[0].getIndex());
        assertEquals(expectedCT[1].getIndex(),observedCT[1].getIndex());

        PhasedGenotypeState GA = new PhasedGenotypeState("GA", "8", 8);
        NucleotideState[] observedGA = getNucleotideState(GA);
        NucleotideState nucleotideStateGA1 = constructor.newInstance("Guanine", "G", 2, (byte)0b0100);
        NucleotideState nucleotideStateGA2 = constructor.newInstance("Adenine", "A", 0, (byte)0b0001);
        NucleotideState[] expectedGA = {nucleotideStateGA1, nucleotideStateGA2};
        assertEquals(expectedGA[0].getIndex(),observedGA[0].getIndex());
        assertEquals(expectedGA[1].getIndex(),observedGA[1].getIndex());

        PhasedGenotypeState GC = new PhasedGenotypeState("GC", "9", 9);
        NucleotideState[] observedGC = getNucleotideState(GC);
        NucleotideState nucleotideStateGC1 = constructor.newInstance("Guanine", "G", 2, (byte)0b0100);
        NucleotideState nucleotideStateGC2 = constructor.newInstance("Cytosine", "C", 1, (byte)0b0010);
        NucleotideState[] expectedGC = {nucleotideStateGC1, nucleotideStateGC2};
        assertEquals(expectedGC[0].getIndex(),observedGC[0].getIndex());
        assertEquals(expectedGC[1].getIndex(),observedGC[1].getIndex());

        PhasedGenotypeState GG = new PhasedGenotypeState("GG", "10", 10);
        NucleotideState[] observedGG = getNucleotideState(GG);
        NucleotideState nucleotideStateGG1 = constructor.newInstance("Guanine", "G", 2, (byte)0b0100);
        NucleotideState nucleotideStateGG2 = constructor.newInstance("Guanine", "G", 2, (byte)0b0100);
        NucleotideState[] expectedGG = {nucleotideStateGG1, nucleotideStateGG2};
        assertEquals(expectedGG[0].getIndex(),observedGG[0].getIndex());
        assertEquals(expectedGG[1].getIndex(),observedGG[1].getIndex());


        PhasedGenotypeState GT = new PhasedGenotypeState("GT", "11", 11);
        NucleotideState[] observedGT = getNucleotideState(GT);
        NucleotideState nucleotideStateGT1 = constructor.newInstance("Guanine", "G", 2, (byte)0b0100);
        NucleotideState nucleotideStateGT2 = constructor.newInstance("Thymine", "T", 3, (byte)0b1000);
        NucleotideState[] expectedGT = {nucleotideStateGT1, nucleotideStateGT2};
        assertEquals(expectedGT[0].getIndex(),observedGT[0].getIndex());
        assertEquals(expectedGT[1].getIndex(),observedGT[1].getIndex());


        PhasedGenotypeState TA = new PhasedGenotypeState("TA", "12", 12);
        NucleotideState[] observedTA = getNucleotideState(TA);
        NucleotideState nucleotideStateTA1 = constructor.newInstance("Thymine", "T", 3, (byte)0b1000);
        NucleotideState nucleotideStateTA2 = constructor.newInstance("Adenine", "A", 0, (byte)0b0001);
        NucleotideState[] expectedTA = {nucleotideStateTA1, nucleotideStateTA2};
        assertEquals(expectedTA[0].getIndex(),observedTA[0].getIndex());
        assertEquals(expectedTA[1].getIndex(),observedTA[1].getIndex());

        PhasedGenotypeState TC = new PhasedGenotypeState("TC", "13", 13);
        NucleotideState[] observedTC = getNucleotideState(TC);
        NucleotideState nucleotideStateTC1 = constructor.newInstance("Thymine", "T", 3, (byte)0b1000);
        NucleotideState nucleotideStateTC2 = constructor.newInstance("Cytosine", "C", 1, (byte)0b0010);
        NucleotideState[] expectedTC = {nucleotideStateTC1, nucleotideStateTC2};
        assertEquals(expectedTC[0].getIndex(),observedTC[0].getIndex());
        assertEquals(expectedTC[1].getIndex(),observedTC[1].getIndex());

        PhasedGenotypeState TG = new PhasedGenotypeState("TG", "14", 14);
        NucleotideState[] observedTG = getNucleotideState(TG);
        NucleotideState nucleotideStateTG1 = constructor.newInstance("Thymine", "T", 3, (byte)0b1000);
        NucleotideState nucleotideStateTG2 = constructor.newInstance("Guanine", "G", 2, (byte)0b0100);
        NucleotideState[] expectedTG = {nucleotideStateTG1, nucleotideStateTG2};
        assertEquals(expectedTG[0].getIndex(),observedTG[0].getIndex());
        assertEquals(expectedTG[1].getIndex(),observedTG[1].getIndex());

        PhasedGenotypeState TT = new PhasedGenotypeState("TT", "15", 15);
        NucleotideState[] observedTT = getNucleotideState(TT);
        NucleotideState nucleotideStateTT1 = constructor.newInstance("Thymine", "T", 3, (byte)0b1000);
        NucleotideState nucleotideStateTT2 = constructor.newInstance("Thymine", "T", 3, (byte)0b1000);
        NucleotideState[] expectedTT = {nucleotideStateTT1, nucleotideStateTT2};
        assertEquals(expectedTT[0].getIndex(),observedTT[0].getIndex());
        assertEquals(expectedTT[1].getIndex(),observedTT[1].getIndex());

        PhasedGenotypeState M = new PhasedGenotypeState("ac", "M", 16);
        NucleotideState[] observedM = getNucleotideState(M);
        NucleotideState nucleotideStateM1 = constructor.newInstance("A/C", "M", 6, (byte)0b0011);
        NucleotideState nucleotideStateM2 = constructor.newInstance("A/C", "M", 6, (byte)0b0011);
        NucleotideState[] expectedM = {nucleotideStateM1, nucleotideStateM2};
        assertEquals(expectedM[0].getIndex(),observedM[0].getIndex());
        assertEquals(expectedM[1].getIndex(),observedM[1].getIndex());

        PhasedGenotypeState R = new PhasedGenotypeState("ag", "R", 17);
        NucleotideState[] observedR = getNucleotideState(R);
        NucleotideState nucleotideStateR1 = constructor.newInstance("A/G", "R", 4, (byte)0b0101);
        NucleotideState nucleotideStateR2 = constructor.newInstance("A/G", "R", 4, (byte)0b0101);
        NucleotideState[] expectedR = {nucleotideStateR1, nucleotideStateR2};
        assertEquals(expectedR[0].getIndex(),observedR[0].getIndex());
        assertEquals(expectedR[1].getIndex(),observedR[1].getIndex());

        PhasedGenotypeState W = new PhasedGenotypeState("at", "W", 18);
        NucleotideState[] observedW = getNucleotideState(W);
        NucleotideState nucleotideStateW1 = constructor.newInstance("A/T", "W", 7, (byte)0b1001);
        NucleotideState nucleotideStateW2 = constructor.newInstance("A/T", "W", 7, (byte)0b1001);
        NucleotideState[] expectedW = {nucleotideStateW1, nucleotideStateW2};
        assertEquals(expectedW[0].getIndex(),observedW[0].getIndex());
        assertEquals(expectedW[1].getIndex(),observedW[1].getIndex());

        PhasedGenotypeState S = new PhasedGenotypeState("cg", "S", 19);
        NucleotideState[] observedS = getNucleotideState(S);
        NucleotideState nucleotideStateS1 = constructor.newInstance("C/G", "S", 8, (byte)0b0110);
        NucleotideState nucleotideStateS2 = constructor.newInstance("C/G", "S", 8, (byte)0b0110);
        NucleotideState[] expectedS = {nucleotideStateS1, nucleotideStateS2};
        assertEquals(expectedS[0].getIndex(),observedS[0].getIndex());
        assertEquals(expectedS[1].getIndex(),observedS[1].getIndex());

        PhasedGenotypeState Y = new PhasedGenotypeState("ct", "Y", 20);
        NucleotideState[] observedY = getNucleotideState(Y);
        NucleotideState nucleotideStateY1 = constructor.newInstance("C/T", "Y", 5, (byte)0b1010);
        NucleotideState nucleotideStateY2 = constructor.newInstance("C/T", "Y", 5, (byte)0b1010);
        NucleotideState[] expectedY = {nucleotideStateY1, nucleotideStateY2};
        assertEquals(expectedY[0].getIndex(),observedY[0].getIndex());
        assertEquals(expectedY[1].getIndex(),observedY[1].getIndex());

        PhasedGenotypeState K = new PhasedGenotypeState("gt", "K", 21);
        NucleotideState[] observedK = getNucleotideState(K);
        NucleotideState nucleotideStateK1 = constructor.newInstance("G/T", "K", 9, (byte)0b1100);
        NucleotideState nucleotideStateK2 = constructor.newInstance("G/T", "K", 9, (byte)0b1100);
        NucleotideState[] expectedK = {nucleotideStateK1, nucleotideStateK2};
        assertEquals(expectedK[0].getIndex(),observedK[0].getIndex());
        assertEquals(expectedK[1].getIndex(),observedK[1].getIndex());

        PhasedGenotypeState unkown = new PhasedGenotypeState("unknown genotype", "?", 22);
        NucleotideState[] observedUnkown = getNucleotideState(unkown);
        NucleotideState nucleotideStateUnkown1 = constructor.newInstance("Unknown base", "?", 15, (byte)0b1111);
        NucleotideState nucleotideStateUnkown2 = constructor.newInstance("Unknown base", "?", 15, (byte)0b1111);
        NucleotideState[] expectedUnkown = {nucleotideStateUnkown1, nucleotideStateUnkown2};
        assertEquals(expectedUnkown[0].toString(),observedUnkown[0].toString());
        assertEquals(expectedUnkown[1].toString(),observedUnkown[1].toString());

        PhasedGenotypeState gap = new PhasedGenotypeState("gap", "-", 23);
        NucleotideState[] observedGap = getNucleotideState(gap);
        NucleotideState nucleotideStateGap1 = constructor.newInstance("Gap", "-", 16, (byte)0b1111);
        NucleotideState nucleotideStateGap2 = constructor.newInstance("Gap", "-", 16, (byte)0b1111);
        NucleotideState[] expectedGap = {nucleotideStateGap1, nucleotideStateGap2};
        assertEquals(expectedGap[0].toString(),observedGap[0].toString());
        assertEquals(expectedGap[1].toString(),observedGap[1].toString());
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
        NucleotideState parent1_state = constructor.newInstance("name", "code", -1 , (byte) 0b11111);
        NucleotideState parent2_state = constructor.newInstance("name", "code", -1 , (byte) 0b11111) ;
        NucleotideState parentA = constructor.newInstance("Adenine", "A", 0, (byte)0b0001);
        NucleotideState parentC = constructor.newInstance("Cytosine", "C",1, (byte)0b0010);
        NucleotideState parentG = constructor.newInstance("Guanine", "G", 2, (byte)0b0100);
        NucleotideState parentT = constructor.newInstance("Thymine", "T", 3, (byte)0b1000);

        // start to test
        // parent1 = A
        parent1_state = parentA;
        parent2_state = parentA;
        PhasedGenotypeState observedAA = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedAA = new PhasedGenotypeState("AA", "0", 0);
        assertEquals(expectedAA.toString(), observedAA.toString());

        parent2_state = parentC;
        PhasedGenotypeState observedAC = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedAC = new PhasedGenotypeState("AC", "1", 1);
        assertEquals(expectedAC.toString(), observedAC.toString());

        parent2_state = parentG;
        PhasedGenotypeState observedAG = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedAG = new PhasedGenotypeState("AG", "2", 2);
        assertEquals(expectedAG.toString(), observedAG.toString());


        parent2_state = parentT;
        PhasedGenotypeState observedAT = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedAT = new PhasedGenotypeState("AT", "3", 3);
        assertEquals(expectedAT.toString(), observedAT.toString());

        // parent1 = C
        parent1_state = parentC;
        parent2_state = parentA;
        PhasedGenotypeState observedCA = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedCA = new PhasedGenotypeState("CA", "4", 4);
        assertEquals(expectedCA.toString(), observedCA.toString());

        parent2_state = parentC;
        PhasedGenotypeState observedCC = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedCC = new PhasedGenotypeState("CC", "5", 5);
        assertEquals(expectedCC.toString(), observedCC.toString());

        parent2_state = parentG;
        PhasedGenotypeState observedCG = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedCG = new PhasedGenotypeState("CG", "6", 6);
        assertEquals(expectedCG.toString(), observedCG.toString());

        parent2_state = parentT;
        PhasedGenotypeState observedCT = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedCT = new PhasedGenotypeState("CT", "7", 7);
        assertEquals(expectedCT.toString(), observedCT.toString());

        // parent1 = G
        parent1_state = parentG;
        parent2_state = parentA;
        PhasedGenotypeState observedGA = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedGA = new PhasedGenotypeState("GC", "8", 8);
        assertEquals(expectedGA.toString(), observedGA.toString());

        parent2_state = parentC;
        PhasedGenotypeState observedGC = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedGC = new PhasedGenotypeState("GC", "9", 9);
        assertEquals(expectedGC.toString(), observedGC.toString());

        parent2_state = parentG;
        PhasedGenotypeState observedGG = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedGG = new PhasedGenotypeState("GG", "a", 10);
        assertEquals(expectedGG.toString(), observedGG.toString());

        parent2_state = parentT;
        PhasedGenotypeState observedGT = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedGT = new PhasedGenotypeState("GT", "b", 11);
        assertEquals(expectedGT.toString(), observedGT.toString());

        // parent1 = T
        parent1_state = parentT;
        parent2_state = parentA;
        PhasedGenotypeState observedTA = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedTA = new PhasedGenotypeState("TA", "c", 12);
        assertEquals(expectedTA.toString(), observedTA.toString());

        parent2_state = parentC;
        PhasedGenotypeState observedTC = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedTC = new PhasedGenotypeState("TC", "d", 13);
        assertEquals(expectedTC.toString(), observedTC.toString());

        parent2_state = parentG;
        PhasedGenotypeState observedTG = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedTG = new PhasedGenotypeState("TG", "e", 14);
        assertEquals(expectedTG.toString(), observedTG.toString());

        parent2_state = parentT;
        PhasedGenotypeState observedTT = getPhasedGenotypeState(parent1_state, parent2_state);
        PhasedGenotypeState expectedTT = new PhasedGenotypeState("TT", "f", 15);
        assertEquals(expectedTT.toString(), observedTT.toString());
    }
}