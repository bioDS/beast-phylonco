package phylonco.beast.evolution.datatype;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;

@Description("Ternary data type for data with three states")
public class Ternary extends DataType.Base  {

    int[][] x = {
            {0},
            {1},
            {2},
            {0, 1, 2}, // gap -
            {0, 1, 2} // missing ?
    };

    public Ternary() {
        stateCount = 3;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "012" + GAP_CHAR + MISSING_CHAR;
    }

    @Override
    public String getTypeDescription() {
        return "ternary";
    }
}
