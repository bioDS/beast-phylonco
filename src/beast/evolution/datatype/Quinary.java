package beast.evolution.datatype;

public class Quinary extends DataType.Base {

    int[][] x = {
            {0},
            {1},
            {2},
            {3},
            {4},
            {0, 1, 2, 3, 4},
    };

    public Quinary() {
        stateCount = 5;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "012345" + MISSING_CHAR;
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();
    }

    @Override
    public String getTypeDescription() {
        return "quinary";
    }

}
