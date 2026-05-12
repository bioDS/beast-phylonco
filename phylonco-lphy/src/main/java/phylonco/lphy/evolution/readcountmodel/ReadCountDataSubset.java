package phylonco.lphy.evolution.readcountmodel;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.Arrays;

public class ReadCountDataSubset extends DeterministicFunction<ReadCountData> {
    private int startSite;
    private int endSite;
    private Value<ReadCountData> readCountData;


    public ReadCountDataSubset(
            @ParameterInfo(name = "rc", narrativeName = "read count data", description = "input data read count.") Value<ReadCountData> readCountData,
            @ParameterInfo(name = "start", narrativeName = "start site index", description = "the start site index to extra data") Value<Integer> startSite,
            @ParameterInfo(name = "end", narrativeName = "end site index", description = "the start site index to extra data") Value<Integer> endSite) {
        if (readCountData == null) {
            throw new IllegalArgumentException("readCount cannot be null");
        }
        if (startSite == null) {
            throw new IllegalArgumentException("start site cannot be null");
        }
        if (endSite == null) {
            throw new IllegalArgumentException("end site cannot be null");
        }
        this.readCountData = readCountData;
        this.startSite = startSite.value();
        this.endSite = endSite.value();
    }

    @GeneratorInfo(name="readCountDataSubset",
            narrativeName = "read count data subset",
            description = "get a read count data subset by 2 site indexes.")

    @Override
    public Value<ReadCountData> apply() {
        ReadCountData rc = readCountData.value();
        String[] chrom = Arrays.copyOfRange(rc.getChromNames(), startSite, endSite+1);
        int[] ref = Arrays.copyOfRange(rc.getRefIndex(), startSite, endSite+1);
        int[] site = Arrays.copyOfRange(rc.getSitesIndex(), startSite, endSite+1);
        ReadCount[][] readCountDataMatrix = new ReadCount[rc.getTaxa().getDimension()][ref.length];
        for (int i = 0; i < rc.getTaxa().getDimension(); i++) {
            for (int j = 0; j < ref.length; j++) {
                readCountDataMatrix[i][j] = rc.getReadCountDataMatrix()[i][startSite+j];
            }
        }
        ReadCountData subsetReadCountData = new ReadCountData(chrom, ref, rc.getTaxa(), readCountDataMatrix, site);
        return new Value<>(null, subsetReadCountData, this);
    }
}
