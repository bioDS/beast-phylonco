package phylonco.lphy.evolution.alignment;

import org.junit.Test;

import static org.junit.Assert.assertArrayEquals;
import static org.junit.Assert.assertEquals;
import static phylonco.lphy.evolution.datatype.NucleotideGenotypeHelper.getNucleotideIndex;
import static phylonco.lphy.evolution.datatype.NucleotideGenotypeHelper.getPhasedGenotypeIndex;

public class NucleotideGenotypeHelperTest {
    @Test
    public void getNucleotideIndexTest() {
        int stateIndex = 0;
        for (int i = 0; i<24;i++){
            int[] observed = getNucleotideIndex(stateIndex);
            if (stateIndex == 0){
                int[] expected = {0,0};
                assertArrayEquals("StateIndex is 0",expected,observed);
            } else if (stateIndex == 1){
                int[] expected = {0,1};
                assertArrayEquals("StateIndex is 1",expected,observed);
            } else if (stateIndex == 2){
                int[] expected = {0,2};
                assertArrayEquals("StateIndex is 2",expected,observed);
            } else if (stateIndex == 3){
                int[] expected = {0,3};
                assertArrayEquals("StateIndex is 3",expected,observed);
            } else if (stateIndex == 4){
                int[] expected = {1,0};
                assertArrayEquals("StateIndex is 4",expected,observed);
            } else if (stateIndex == 5){
                int[] expected = {1,1};
                assertArrayEquals("StateIndex is 5",expected,observed);
            } else if (stateIndex == 6){
                int[] expected = {1,2};
                assertArrayEquals("StateIndex is 6",expected,observed);
            } else if (stateIndex == 7){
                int[] expected = {1,3};
                assertArrayEquals("StateIndex is 7",expected,observed);
            } else if (stateIndex == 8){
                int[] expected = {2,0};
                assertArrayEquals("StateIndex is 8",expected,observed);
            } else if (stateIndex == 9){
                int[] expected = {2,1};
                assertArrayEquals("StateIndex is 9",expected,observed);
            } else if (stateIndex == 10){
                int[] expected = {2,2};
                assertArrayEquals("StateIndex is 10",expected,observed);
            } else if (stateIndex == 11){
                int[] expected = {2,3};
                assertArrayEquals("StateIndex is 11",expected,observed);
            } else if (stateIndex == 12){
                int[] expected = {3,0};
                assertArrayEquals("StateIndex is 12",expected,observed);
            } else if (stateIndex == 13){
                int[] expected = {3,1};
                assertArrayEquals("StateIndex is 13",expected,observed);
            } else if (stateIndex == 14){
                int[] expected = {3,2};
                assertArrayEquals("StateIndex is 14",expected,observed);
            } else if (stateIndex == 15){
                int[] expected = {3,3};
                assertArrayEquals("StateIndex is 15",expected,observed);
            } else if (stateIndex == 16){
                int[] expected = {6,6};
                assertArrayEquals("StateIndex is 16",expected,observed);
            } else if (stateIndex == 17){
                int[] expected = {4,4};
                assertArrayEquals("StateIndex is 17",expected,observed);
            } else if (stateIndex == 18){
                int[] expected = {7,7};
                assertArrayEquals("StateIndex is 17",expected,observed);
            } else if (stateIndex == 19){
                int[] expected = {8,8};
                assertArrayEquals("StateIndex is 19",expected,observed);
            } else if (stateIndex == 20){
                int[] expected = {5,5};
                assertArrayEquals("StateIndex is 20",expected,observed);
            } else if (stateIndex == 21){
                int[] expected = {9,9};
                assertArrayEquals("StateIndex is 21",expected,observed);
            } else if (stateIndex == 22){
                int[] expected = {15,15};
                assertArrayEquals("StateIndex is 22",expected,observed);
            } else if (stateIndex == 23){
                int[] expected = {16,16};
                assertArrayEquals("StateIndex is 23",expected,observed);
            }
            stateIndex++;
        }
    }

    @Test
    public void getPhasedGenotypeIndexTest() {
        int parent1_index;
        int parent2_index;
        for ( parent1_index = 0; parent1_index < 4;parent1_index++){
            for ( parent2_index = 0; parent2_index < 4; parent2_index++){
                int observed = getPhasedGenotypeIndex(parent1_index,parent2_index);

                if (parent1_index == 0){
                    if (parent2_index == 0){
                        int expected = 0;
                        assertEquals("0-0",expected,observed);
                    } else if (parent2_index == 1){
                        int expected = 1;
                        assertEquals("0-1",expected,observed);
                    } else if (parent2_index == 2){
                        int expected = 2;
                        assertEquals("0-2",expected,observed);
                    } else if (parent2_index == 3){
                        int expected = 3;
                        assertEquals("0-3",expected,observed);
                    }
                } else if (parent1_index == 1){
                    if (parent2_index == 0){
                        int expected = 4;
                        assertEquals("1-0",expected,observed);
                    } else if (parent2_index == 1){
                        int expected = 5;
                        assertEquals("1-1",expected,observed);
                    } else if (parent2_index == 2){
                        int expected = 6;
                        assertEquals("1-2",expected,observed);
                    } else if (parent2_index == 3){
                        int expected = 7;
                        assertEquals("1-3",expected,observed);
                    }
                } else if (parent1_index == 2){
                    if (parent2_index == 0){
                        int expected = 8;
                        assertEquals("2-0",expected,observed);
                    } else if (parent2_index == 1){
                        int expected = 9;
                        assertEquals("2-1",expected,observed);
                    } else if (parent2_index == 2){
                        int expected = 10;
                        assertEquals("2-2",expected,observed);
                    } else if (parent2_index == 3){
                        int expected = 11;
                        assertEquals("2-3",expected,observed);
                    }
                } else if (parent1_index == 3){
                    if (parent2_index == 0){
                        int expected = 12;
                        assertEquals("3-0",expected,observed);
                    } else if (parent2_index == 1){
                        int expected = 13;
                        assertEquals("3-1",expected,observed);
                    } else if (parent2_index == 2){
                        int expected = 14;
                        assertEquals("3-2",expected,observed);
                    } else if (parent2_index == 3){
                        int expected = 15;
                        assertEquals("3-3",expected,observed);
                    }
                }
            }
        }
    }
}
