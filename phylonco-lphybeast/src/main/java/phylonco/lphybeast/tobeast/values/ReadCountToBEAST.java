package phylonco.lphybeast.tobeast.values;

import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import phylonco.beast.evolution.datatype.ReadCount;
import phylonco.lphy.evolution.readcountmodel.ReadCountData;

public class ReadCountToBEAST implements ValueToBEAST <ReadCountData, ReadCount> {

    @Override
    public ReadCount valueToBEAST(Value<ReadCountData> value, BEASTContext context) {
//        ReadCount readCount = new ReadCount(value.value().getTaxa().getDimension(), value.value().nchar());
        ReadCount readCount = new ReadCount();
        String readC = "\n";
        int n = value.value().getTaxa().getDimension();
        int l = value.value().nchar();
        for (int i = 0; i < n; i++) {
            // n taxa
            for (int j = 0; j < l; j++) {
                // n sites
                int countA = value.value().getState(i,j).getCount("A");
                int countC = value.value().getState(i,j).getCount("C");
                int countG = value.value().getState(i,j).getCount("G");
                int countT = value.value().getState(i,j).getCount("T");
                readC += String.format("%d,%d,%d,%d; ", countA, countC, countG, countT);
            }
            readC += "\n";
        }

//        for (int i = 0; i < value.value().getTaxa().getDimension(); i++) {
//            for (int j = 0; j < value.value().nchar(); j++) {
//                readCount.setReadCounts(i,j,value.value().getState(i,j).getReadCounts());
//            }
//        }
        readCount.setInputValue("value", readC);
        readCount.initAndValidate();
        return readCount;
    }

    @Override
    public Class getValueClass() {
        return ReadCountData.class;
    }
}