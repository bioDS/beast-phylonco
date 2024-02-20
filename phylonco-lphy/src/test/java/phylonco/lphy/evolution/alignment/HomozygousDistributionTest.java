package phylonco.lphy.evolution.alignment;

import jebl.evolution.align.Align;
import jebl.evolution.sequences.SequenceType;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.core.model.Value;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class HomozygousDistributionTest {
    /**
     * Test getStateIndex function in {@link HomozygousAlignmentDistribution}
     * @Parameters the nucleotide state index
     * @Return the homozygous state phased genotype index of the input nucleotide state
     */
    @Test
    public void getAmbiguousStateIndexTest() {
        Alignment alignment = new SimpleAlignment(Taxa.createTaxa(10),10, SequenceType.NUCLEOTIDE);
        Value<Alignment> alignmentValue = new Value<>("alignment", alignment);
        int originalStateIndex;
        HomozygousAlignmentDistribution instance = new HomozygousAlignmentDistribution(alignmentValue);
        for (originalStateIndex = 4; originalStateIndex < 17; originalStateIndex++){
            if (originalStateIndex == 4){
                int R_A = 0;
                int R_G = 0;
                int R_exc = 0;
                for (int i = 0; i < 10000; i++){
                    int certainState_R = instance.getAmbiguousStateIndex(originalStateIndex);
                    if (certainState_R == 0){
                        R_A++;
                    } else if (certainState_R == 2){
                        R_G++;
                    } else {
                        R_exc++;
                    }
                }
                assertEquals("For R_A",5000 , R_A , 500);
                assertEquals("For R_G",5000 , R_G , 500);
                assertEquals("For R Exceptions",0 , R_exc , 0);
            } else if (originalStateIndex == 5){
                int Y_C = 0;
                int Y_T = 0;
                int Y_exc = 0;
                for (int i = 0; i < 10000; i++){
                    int certainState_Y = instance.getAmbiguousStateIndex(originalStateIndex);
                    if (certainState_Y == 1){
                        Y_C++;
                    } else if (certainState_Y == 3){
                        Y_T++;
                    } else {
                        Y_exc++;
                    }
                }
                assertEquals("For Y_C",5000 , Y_C , 500);
                assertEquals("For Y_T",5000 , Y_T , 500);
                assertEquals("For Y Exceptions",0 , Y_exc , 0);
            } else if (originalStateIndex == 6){
                int M_A = 0;
                int M_C = 0;
                int M_exc = 0;
                for (int i = 0; i < 10000; i++){
                    int certainState_M = instance.getAmbiguousStateIndex(originalStateIndex);

                    if (certainState_M == 0){
                        M_A++;
                    } else if (certainState_M == 1){
                        M_C++;
                    } else {
                        M_exc++;
                    }
                }

                assertEquals("For M_A",5000 , M_A , 500);
                assertEquals("For M_C",5000 , M_C , 500);
                assertEquals("For M Exceptions",0 , M_exc , 0);
            } else if (originalStateIndex == 7){
                int W_A = 0;
                int W_T = 0;
                int W_exc = 0;
                for (int i = 0; i < 10000; i++){
                    int certainState_W = instance.getAmbiguousStateIndex(originalStateIndex);

                    if (certainState_W == 0){
                        W_A++;
                    } else if (certainState_W == 3){
                        W_T++;
                    } else {
                        W_exc++;
                    }
                }
                assertEquals("For W_A",5000 , W_A , 500);
                assertEquals("For W_T",5000 , W_T , 500);
                assertEquals("For W Exceptions",0 , W_exc , 0);
            } else if (originalStateIndex == 8){
                int S_C = 0;
                int S_G = 0;
                int S_exc = 0;
                for (int i = 0; i < 10000; i++){
                    int certainState_S = instance.getAmbiguousStateIndex(originalStateIndex);

                    if (certainState_S == 1){
                        S_C++;
                    } else if (certainState_S == 2){
                        S_G++;
                    } else {
                        S_exc++;
                    }
                }
                assertEquals("For S_C",5000 , S_C , 500);
                assertEquals("For S_G", 5000 , S_G , 500);
                assertEquals("For S Exceptions",0 , S_exc , 0);
            } else if (originalStateIndex == 9){
                int K_G = 0;
                int K_T = 0;
                int K_exc = 0;
                for (int i = 0; i < 10000; i++){
                    int certainState_K = instance.getAmbiguousStateIndex(originalStateIndex);

                    if (certainState_K == 2){
                        K_G++;
                    } else if (certainState_K == 3){
                        K_T++;
                    } else {
                        K_exc++;
                    }
                }
                assertEquals("For K_G",5000 , K_G , 500);
                assertEquals("For K_T", 5000, K_T , 500);
                assertEquals("For K Exceptions",0 , K_exc , 0);
            } else if (originalStateIndex == 10){
                int B_C = 0;
                int B_G = 0;
                int B_T = 0;
                int B_exc = 0;
                for (int i = 0; i < 10000; i++){
                    int certainState_B = instance.getAmbiguousStateIndex(originalStateIndex);

                    if (certainState_B == 1){
                        B_C++;
                    } else if (certainState_B == 2){
                        B_G++;
                    } else if (certainState_B == 3){
                        B_T++;
                    }else {
                        B_exc++;
                    }
                }
                assertEquals("For B_C",3333 , B_C , 333);
                assertEquals("For B_G",3333 , B_G , 333);
                assertEquals("For B_T",3333 , B_T , 333);
                assertEquals("For B Exceptions",0 , B_exc , 0);
            } else if (originalStateIndex == 11){
                int D_A = 0;
                int D_G = 0;
                int D_T = 0;
                int D_exc = 0;
                for (int i = 0; i < 10000; i++){
                    int certainState_D = instance.getAmbiguousStateIndex(originalStateIndex);

                    if (certainState_D == 0){
                        D_A++;
                    } else if (certainState_D == 2){
                        D_G++;
                    } else if (certainState_D == 3){
                        D_T++;
                    }else {
                        D_exc++;
                    }
                }
                assertEquals("For D_A",3333 , D_A , 333);
                assertEquals("For D_G",3333 , D_G , 333);
                assertEquals("For D_T",3333 , D_T , 333);
                assertEquals("For D Exceptions",0 , D_exc , 0);
            } else if (originalStateIndex == 12){
                int H_A = 0;
                int H_C = 0;
                int H_T = 0;
                int H_exc = 0;
                for (int i = 0; i < 10000; i++){
                    int certainState_H = instance.getAmbiguousStateIndex(originalStateIndex);

                    if (certainState_H == 0){
                        H_A++;
                    } else if (certainState_H == 1){
                        H_C++;
                    } else if (certainState_H == 3){
                        H_T++;
                    }else {
                        H_exc++;
                    }
                }
                assertEquals("For H_A",3333 , H_A , 333);
                assertEquals("For H_C",3333 , H_C , 333);
                assertEquals("For H_T",3333 , H_T , 333);
                assertEquals("For H Exceptions",0 , H_exc , 0);
            } else if (originalStateIndex == 13){
                int V_A = 0;
                int V_C = 0;
                int V_G = 0;
                int V_exc = 0;
                for (int i = 0; i < 10000; i++){
                    int certainState_V = instance.getAmbiguousStateIndex(originalStateIndex);

                    if (certainState_V == 0){
                        V_A++;
                    } else if (certainState_V == 1){
                        V_C++;
                    } else if (certainState_V == 2){
                        V_G++;
                    }else {
                        V_exc++;
                    }
                }
                assertEquals("For V_A",3333 , V_A , 333);
                assertEquals("For V_C",3333 , V_C , 333);
                assertEquals("For V_G",3333 , V_G , 333);
                assertEquals("For V Exceptions",0 , V_exc , 0);
            } else if (originalStateIndex >=14 && originalStateIndex<17){
                int num_A = 0;
                int num_C = 0;
                int num_G = 0;
                int num_T = 0;
                int num_exc = 0;
                for (int i = 0; i < 10000; i++){
                    int certainState = instance.getAmbiguousStateIndex(originalStateIndex);
                    if (certainState == 0){
                        num_A++;
                    } else if (certainState == 1){
                        num_C++;
                    } else if (certainState == 2){
                        num_G++;
                    }else if (certainState == 3){
                        num_T++;
                    } else {
                        num_exc++;
                    }
                }
                assertEquals("For Unkown States A",2500 , num_A , 250);
                assertEquals("For Unkown States C",2500 , num_C , 250);
                assertEquals("For Unkown States G",2500 , num_G , 250);
                assertEquals("For Unkown States T",2500 , num_T , 250);
                assertEquals("For Unkown States Exceptions",0 , num_exc , 0);
            }
        }
    }
}
