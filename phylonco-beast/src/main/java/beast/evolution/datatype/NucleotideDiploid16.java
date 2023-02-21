package beast.evolution.datatype;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;

@Description("Phased diploid nucleotide data type")
public class NucleotideDiploid16 extends DataType.Base {
    int[][] x = {
            {0}, // AA - 0
            {1}, // AC - 1
            {2}, // AG - 2
            {3}, // AT - 3
            {4}, // CA - 4
            {5}, // CC - 5
            {6}, // CG - 6
            {7}, // CT - 7
            {8}, // GA - 8
            {9}, // GC - 9
            {10}, // GG - a
            {11}, // GT - b
            {12}, // TA - c
            {13}, // TC - d
            {14}, // TG - e
            {15}, // TT - f
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
        codeMap = "0123456789ABCDEFMRWSYK" + GAP_CHAR + MISSING_CHAR;
    }

    @Override
    public String getTypeDescription() {
        return "nucleotideDiploid16";
    }

}