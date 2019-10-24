package beast.evolution.datatype;

import beast.evolution.datatype.DataType;

public class NucleotideDiploid extends DataType.Base {
    int[][] x = {
            {0}, // AA - A
            {1}, // AC - M
            {2}, // AG - R
            {3}, // AT - W
            {4}, // CC - C
            {5}, // CG - S
            {6}, // CT - Y
            {7}, // GG - G
            {8}, // GT - K
            {9}, // TT - T
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9} // any base - ?
    };

    public NucleotideDiploid() {
        stateCount = 10;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "AMRWCSYGKT" + MISSING_CHAR;
    }

    @Override
    public String getTypeDescription() {
        return "nucleotideDiploid";
    }

}
