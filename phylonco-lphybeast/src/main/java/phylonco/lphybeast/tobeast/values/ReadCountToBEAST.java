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
        String readC = new String();
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
                readC += String.format("%d:%d:%d:%d", countA, countC, countG, countT);
                if (j < l - 1) {
                    readC += ",";
                }
                // read count per site
            }
            readC += "\n"; // new line character per taxa
        }

        String[] taxaNames = value.value().getTaxaNames();
        String taxaString = new String();
        for (int i = 0; i < taxaNames.length; i++) {
            taxaString += taxaNames[i];
            taxaString += " ";
        }



        readCount.setInputValue("value", readC);
        readCount.setInputValue("taxaNames", taxaString);
        readCount.initAndValidate();
        return readCount;
    }

    @Override
    public Class getValueClass() {
        return ReadCountData.class;
    }
}
