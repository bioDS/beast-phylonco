package datatype;

import beast.evolution.datatype.DataType;

public class TernaryWithError extends DataType.Base implements DataTypeWithError {

    int[][] x = {
            {0},  // 0 Homozygous reference
            {1},  // 1 Heterozygous
            {2},  // 2 Homozygous non reference
            {0, 1, 2}, // ?
    };

    public TernaryWithError() {
        stateCount = 3;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "012" + MISSING_CHAR;
    }

    @Override
    public String getTypeDescription() {
        return "ternaryWithError";
    }
}
