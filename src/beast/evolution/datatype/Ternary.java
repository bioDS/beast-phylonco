package beast.evolution.datatype;

public class Ternary extends DataType.Base  {

    int[][] x = {
            {0},
            {1},
            {2},
            {0, 1, 2},
    };

    public Ternary() {
        stateCount = 3;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "012" + MISSING_CHAR;
    }

    @Override
    public String getTypeDescription() {
        return "ternary";
    }
}
