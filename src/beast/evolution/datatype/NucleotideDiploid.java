package beast.evolution.datatype;

import beast.evolution.datatype.DataType;

public class NucleotideDiploid extends DataType.Base {
    int[][] x = {
            {0}, // AA - A
            {1}, // CC - C
            {2}, // GG - G
            {3}, // TT - T
            {4}, // AT - W
            {5}, // CG - S
            {6}, // AC - M
            {7}, // GT - K
            {8}, // AG - R
            {9}, // CT - Y
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9} // any base - ?
    };

    public NucleotideDiploid() {
        stateCount = 10;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "ACGTWSMKRY" + MISSING_CHAR;
    }

    @Override
    public String getTypeDescription() {
        return "nucleotideDiploid";
    }

}
