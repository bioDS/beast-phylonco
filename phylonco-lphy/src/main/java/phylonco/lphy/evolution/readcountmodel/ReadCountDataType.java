package phylonco.lphy.evolution.readcountmodel;

import lphy.core.model.Value;

public class ReadCountDataType extends Value<ReadCountData> {

    public ReadCountDataType(String id, ReadCountData value) {
        super(id, value);
    }

    @Override
    public String toString() {
        String result = "";
        int n = value.getTaxa().getDimension();;
        int l = value.nchar();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < l; j++) {
                int countA = value.readCountDataMatrix[i][j].getCount("A");
                int countC = value.readCountDataMatrix[i][j].getCount("C");
                int countG = value.readCountDataMatrix[i][j].getCount("G");
                int countT = value.readCountDataMatrix[i][j].getCount("T");
                result += String.format("A: %d, C: %d, G: %d, T: %d; \t", countA, countC, countG, countT);
            }
        }
        return result;
    }
}
