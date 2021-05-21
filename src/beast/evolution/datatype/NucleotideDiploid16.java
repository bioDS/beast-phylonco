package beast.evolution.datatype;

import beast.core.Description;

@Description("Phased diploid nucleotide data type")
public class NucleotideDiploid16 extends DataType.Base {
    int[][] x = {
            {0}, // AA - A
            {1}, // AC - a
            {2}, // AG - b
            {3}, // AT - c
            {4}, // CA - 1
            {5}, // CC - C
            {6}, // CG - d
            {7}, // CT - e
            {8}, // GA - 2
            {9}, // GC - 4
            {10}, // GG - G
            {11}, // GT - f
            {12}, // TA - 3
            {13}, // TC - 5
            {14}, // TG - 6
            {15}, // TT - T
            {1, 4}, // AC or CA - M
            {2, 8}, // AG or GA - R
            {3, 12}, // AT or TA - W
            {6, 9}, // CG or GC - S
            {7, 13}, // CT or TC - Y
            {11, 14}, // GT or TG - K
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, // gap -
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15} // missing ?
    };

    public NucleotideDiploid16() {
        stateCount = 16;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "Aabc1Cde24Gf356TMRWSYK" + GAP_CHAR + MISSING_CHAR;
    }

    @Override
    public String getTypeDescription() {
        return "nucleotideDiploid16";
    }

}