package beast.evolution.datatype;

import beast.core.Description;

@Description("Data type for combined nubleotide and methylation data")
public class NucleotideMethylation extends DataType.Base {
    // TODO currently without all ambiguity codes
    int[][] x = {
            {0},  // A
            {1},  // C
            {2},  // G
            {3},  // T

            {1, 2}, // 0 = C and G
            {4, 5}, // 1 = P and J = C' and G'
            {0, 3}, // W = A and T
            {4},  // P = MetC = C'
            {5},  // J = MetC on opposite strand = G'

            {0, 1, 2, 3, 4, 5}, // N
            {0, 1, 2, 3, 4, 5}, // X
            {0, 1, 2, 3, 4, 5}, // -
            {0, 1, 2, 3, 4, 5}, // ?

            /*
            {3},  // U
            {0, 2}, // R
            {1, 3}, // Y
            {0, 1}, // M
            {0, 3}, // W
            {1, 2}, // S
            {2, 3}, // K
            {1, 2, 3}, // B
            {0, 2, 3}, // D
            {0, 1, 3}, // H
            {0, 1, 2}, // V
            {0, 1, 2, 3}, // N
            {0, 1, 2, 3}, // X
            {0, 1, 2, 3}, // -
            {0, 1, 2, 3}, // ?
            */
    };

    public NucleotideMethylation() {
        stateCount = 6;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "ACGTPJO1WNX" + GAP_CHAR + MISSING_CHAR;
    }

    @Override
    public String getTypeDescription() {
        return "NucleotideMethylation";
    }


    @Override
    public char getChar(int state){
        return codeMap.charAt(state);
    }
}
