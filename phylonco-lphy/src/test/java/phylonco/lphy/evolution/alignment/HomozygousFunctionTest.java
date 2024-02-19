package phylonco.lphy.evolution.alignment;

import lphy.base.distribution.UniformDiscrete;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class HomozygousFunctionTest {
    @Test
    public void getExpectedCanonicalStates() {
        // canonical states --> homozygous states
        int stateIndex_A = 0;
        assertEquals("ForA",0,4 * stateIndex_A + stateIndex_A,0);

        int stateIndex_C = 1;
        assertEquals("ForC",5,4 * stateIndex_C + stateIndex_C,0);

        int stateIndex_G = 2;
        assertEquals("ForG",10,4 * stateIndex_G + stateIndex_G,0);

        int stateIndex_T = 3;
        assertEquals("ForT", 15, 4 * stateIndex_T + stateIndex_T,0);
    }

    @Test
    public void getExpectedR() {
        // ambiguities for state R
        int stateIndex_R = 4;
        int R_A = 0;
        int R_G = 0;
        int R_exc = 0;
        for (int i = 0; i < 10000; i++){
            // get the array for the states
            int[] ambiguousState = ambiguousState(stateIndex_R);
            // get the Value<Integer> for the lower and upper boundary
            Value<Integer> lower = new Value<>("id", 0);
            Value<Integer> upper = new Value<>("id",ambiguousState.length-1);

            // get the random index for the integer in the array
            UniformDiscrete uniformDiscrete = new UniformDiscrete(lower, upper);
            RandomVariable<Integer> randomNumber = uniformDiscrete.sample();

            // give the stateIndex its certain state
            int certainState_R = ambiguousState[randomNumber.value()];

            if (certainState_R == 0){
                R_A++;
            } else if (certainState_R == 2){
                R_G++;
            } else {
                R_exc++;
            }
        }
        assertEquals("For R_A",5000 , R_A , 500);
        assertEquals("For R_G",5000 , R_A , 500);
        assertEquals("For R Exceptions",0 , R_exc , 0);
    }

    @Test
    public void getExpectedY() {
        // ambiguities for state Y
        int stateIndex_Y = 5;
        int Y_C = 0;
        int Y_T = 0;
        int Y_exc = 0;
        for (int i = 0; i < 10000; i++){
            // get the array for the states
            int[] ambiguousState = ambiguousState(stateIndex_Y);
            // get the Value<Integer> for the lower and upper boundary
            Value<Integer> lower = new Value<>("id", 0);
            Value<Integer> upper = new Value<>("id",ambiguousState.length-1);

            // get the random index for the integer in the array
            UniformDiscrete uniformDiscrete = new UniformDiscrete(lower, upper);
            RandomVariable<Integer> randomNumber = uniformDiscrete.sample();

            // give the stateIndex its certain state
            int certainState_Y = ambiguousState[randomNumber.value()];

            if (certainState_Y == 1){
                Y_C++;
            } else if (certainState_Y == 3){
                Y_T++;
            } else {
                Y_exc++;
            }
        }
        assertEquals("For Y_C",5000 , Y_C , 500);
        assertEquals("For Y_T",5000 , Y_C , 500);
        assertEquals("For Y Exceptions",0 , Y_exc , 0);
    }
    @Test
    public void getExpectedM() {
        // ambiguities for state M
        int stateIndex_M = 6;
        int M_A = 0;
        int M_C = 0;
        int M_exc = 0;
        for (int i = 0; i < 10000; i++){
            // get the array for the states
            int[] ambiguousState = ambiguousState(stateIndex_M);
            // get the Value<Integer> for the lower and upper boundary
            Value<Integer> lower = new Value<>("id", 0);
            Value<Integer> upper = new Value<>("id",ambiguousState.length-1);

            // get the random index for the integer in the array
            UniformDiscrete uniformDiscrete = new UniformDiscrete(lower, upper);
            RandomVariable<Integer> randomNumber = uniformDiscrete.sample();

            // give the stateIndex its certain state
            int certainState_M = ambiguousState[randomNumber.value()];

            if (certainState_M == 0){
                M_A++;
            } else if (certainState_M == 1){
                M_C++;
            } else {
                M_exc++;
            }
        }

        assertEquals("For M_A",5000 , M_A , 500);
        assertEquals("For M_C",5000 , M_A , 500);
        assertEquals("For M Exceptions",0 , M_exc , 0);
    }

    @Test
    public void getExpectedW() {
        int stateIndex_W=7;
        int W_A = 0;
        int W_T = 0;
        int W_exc = 0;
        for (int i = 0; i < 10000; i++){
            // get the array for the states
            int[] ambiguousState = ambiguousState(stateIndex_W);
            // get the Value<Integer> for the lower and upper boundary
            Value<Integer> lower = new Value<>("id", 0);
            Value<Integer> upper = new Value<>("id",ambiguousState.length-1);

            // get the random index for the integer in the array
            UniformDiscrete uniformDiscrete = new UniformDiscrete(lower, upper);
            RandomVariable<Integer> randomNumber = uniformDiscrete.sample();

            // give the stateIndex its certain state
            int certainState_W = ambiguousState[randomNumber.value()];

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
    }

    @Test
    public void getExpectedS() {
        int stateIndex_S = 8;
        int S_C = 0;
        int S_G = 0;
        int S_exc = 0;
        for (int i = 0; i < 10000; i++){
            // get the array for the states
            int[] ambiguousState = ambiguousState(stateIndex_S);
            // get the Value<Integer> for the lower and upper boundary
            Value<Integer> lower = new Value<>("id", 0);
            Value<Integer> upper = new Value<>("id",ambiguousState.length-1);

            // get the random index for the integer in the array
            UniformDiscrete uniformDiscrete = new UniformDiscrete(lower, upper);
            RandomVariable<Integer> randomNumber = uniformDiscrete.sample();

            // give the stateIndex its certain state
            int certainState_S = ambiguousState[randomNumber.value()];

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
    }

    @Test
    public void getExpectedK() {
        int stateIndex_K = 9;
        int K_G = 0;
        int K_T = 0;
        int K_exc = 0;
        for (int i = 0; i < 10000; i++){
            // get the array for the states
            int[] ambiguousState = ambiguousState(stateIndex_K);
            // get the Value<Integer> for the lower and upper boundary
            Value<Integer> lower = new Value<>("id", 0);
            Value<Integer> upper = new Value<>("id",ambiguousState.length-1);

            // get the random index for the integer in the array
            UniformDiscrete uniformDiscrete = new UniformDiscrete(lower, upper);
            RandomVariable<Integer> randomNumber = uniformDiscrete.sample();

            // give the stateIndex its certain state
            int certainState_K = ambiguousState[randomNumber.value()];

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
    }

    @Test
    public void getExpectedB() {
        int stateIndex_B = 10;
        int B_C = 0;
        int B_G = 0;
        int B_T = 0;
        int B_exc = 0;
        for (int i = 0; i < 10000; i++){
            // get the array for the states
            int[] ambiguousState = ambiguousState(stateIndex_B);
            // get the Value<Integer> for the lower and upper boundary
            Value<Integer> lower = new Value<>("id", 0);
            Value<Integer> upper = new Value<>("id",ambiguousState.length-1);

            // get the random index for the integer in the array
            UniformDiscrete uniformDiscrete = new UniformDiscrete(lower, upper);
            RandomVariable<Integer> randomNumber = uniformDiscrete.sample();

            // give the stateIndex its certain state
            int certainState_B = ambiguousState[randomNumber.value()];

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
    }

    @Test
    public void getExpectedD() {
        int stateIndex_D = 11;
        int D_A = 0;
        int D_G = 0;
        int D_T = 0;
        int D_exc = 0;
        for (int i = 0; i < 10000; i++){
            // get the array for the states
            int[] ambiguousState = ambiguousState(stateIndex_D);
            // get the Value<Integer> for the lower and upper boundary
            Value<Integer> lower = new Value<>("id", 0);
            Value<Integer> upper = new Value<>("id",ambiguousState.length-1);

            // get the random index for the integer in the array
            UniformDiscrete uniformDiscrete = new UniformDiscrete(lower, upper);
            RandomVariable<Integer> randomNumber = uniformDiscrete.sample();

            // give the stateIndex its certain state
            int certainState_D = ambiguousState[randomNumber.value()];

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
    }

    @Test
    public void getExpectedH() {
        int stateIndex_H = 12;
        int H_A = 0;
        int H_C = 0;
        int H_T = 0;
        int H_exc = 0;
        for (int i = 0; i < 10000; i++){
            // get the array for the states
            int[] ambiguousState = ambiguousState(stateIndex_H);
            // get the Value<Integer> for the lower and upper boundary
            Value<Integer> lower = new Value<>("id", 0);
            Value<Integer> upper = new Value<>("id",ambiguousState.length-1);

            // get the random index for the integer in the array
            UniformDiscrete uniformDiscrete = new UniformDiscrete(lower, upper);
            RandomVariable<Integer> randomNumber = uniformDiscrete.sample();

            // give the stateIndex its certain state
            int certainState_H = ambiguousState[randomNumber.value()];

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
    }

    @Test
    public void getExpectedV() {
        int stateIndex_V = 13;
        int V_A = 0;
        int V_C = 0;
        int V_G = 0;
        int V_exc = 0;
        for (int i = 0; i < 10000; i++){
            // get the array for the states
            int[] ambiguousState = ambiguousState(stateIndex_V);
            // get the Value<Integer> for the lower and upper boundary
            Value<Integer> lower = new Value<>("id", 0);
            Value<Integer> upper = new Value<>("id",ambiguousState.length-1);

            // get the random index for the integer in the array
            UniformDiscrete uniformDiscrete = new UniformDiscrete(lower, upper);
            RandomVariable<Integer> randomNumber = uniformDiscrete.sample();

            // give the stateIndex its certain state
            int certainState_V = ambiguousState[randomNumber.value()];

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
    }

    @Test
    public void getExpectedUnkown() {
        for (int n = 14;n<17;n++){
            int stateIndex_unkown = n;
            int num_A = 0;
            int num_C = 0;
            int num_G = 0;
            int num_T = 0;
            int num_exc = 0;
            for (int i = 0; i < 10000; i++){
                // get the array for the states
                int[] ambiguousState = ambiguousState(stateIndex_unkown);
                // get the Value<Integer> for the lower and upper boundary
                Value<Integer> lower = new Value<>("id", 0);
                Value<Integer> upper = new Value<>("id",ambiguousState.length-1);

                // get the random index for the integer in the array
                UniformDiscrete uniformDiscrete = new UniformDiscrete(lower, upper);
                RandomVariable<Integer> randomNumber = uniformDiscrete.sample();

                // give the stateIndex its certain state
                int certainState = ambiguousState[randomNumber.value()];

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
            assertEquals("For Unkown States C",2500 , num_C , 250);
            assertEquals("For Unkown States G",2500 , num_G , 250);
            assertEquals("For Unkown States T",2500 , num_T , 250);
            assertEquals("For Unkown States Exceptions",0 , num_exc , 0);
        }

    }

    private int[] ambiguousState(int stateIndex) {
        // switch the ambiguous states into canonical states (0=A, 1=C, 2=G, 3=T)
        switch (stateIndex) {
            case 4:
                // 4 = A/G
                return new int[]{0, 2};
            case 5:
                // 5 = C/T
                return new int[]{1, 3};
            case 6:
                // 6 = A/C
                return new int[]{0, 1};
            case 7:
                // 7 = A/T
                return new int[]{0, 3};
            case 8:
                // 8 = C/G
                return new int[]{1, 2};
            case 9:
                // 9 = G/T
                return new int[]{2, 3};
            case 10:
                // 10 = C/G/T
                return new int[]{1, 2, 3};
            case 11:
                // 11 = A/G/T
                return new int[]{0, 2, 3};
            case 12:
                // 12 = A/C/T
                return new int[]{0, 1, 3};
            case 13:
                // 13 = A/C/G
                return new int[]{0, 1, 2};
            case 14, 16, 15:
                // 14 = unkown base (N) = A/C/G/T
                // 15 = unkown base (?) = A/C/G/T
                // 16 = gap (-) = A/C/G/T
                return new int[]{0, 1, 2, 3};
            default:
                throw new RuntimeException("Unexpected state: " + stateIndex);
        }
    }
}
