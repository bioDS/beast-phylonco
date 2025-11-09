package phylonco.lphy.evolution.readcountmodel;

import lphy.base.evolution.CellPosition;
import lphy.base.evolution.Taxa;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.*;

public class CopyReadCount extends DeterministicFunction<ReadCountData> {
    public final String sitesName = "sites";
    public final String readCountName = "readCount";

    public CopyReadCount(@ParameterInfo(name = sitesName, description = "the array sites that wants to copy from the read count data ")Value<CellPosition[]> sites,
                         @ParameterInfo(name = readCountName, description = "the read count data that wants to copy from") Value<ReadCountData> readCount) {
        if (readCount == null) {
            throw new IllegalArgumentException("readCount cannot be null");
        }
        if (sites == null) {
            throw new IllegalArgumentException("sites should at least be 1");
        }

        setParam(sitesName, sites);
        setParam(readCountName, readCount);
    }

    @GeneratorInfo(name = "copyReadCount",
            description = "copy given sites from a read count data")
    @Override
    public Value<ReadCountData> apply() {
        // get params from constructor
        ReadCountData readCount = getReadCount().value();
        CellPosition[] cellPositions = getSites().value();

        // get exclusive taxon and sites
        List<String> taxaNames = new ArrayList<>();
        List<Integer> positions = new ArrayList<>();
        for (int i = 0; i < cellPositions.length; i++) {
            String cellName = cellPositions[i].getCellName();
            int pos = cellPositions[i].getPosition();

            if (!taxaNames.contains(cellName)) {
                taxaNames.add(cellName);
            }
            if (!positions.contains(pos)) {
                positions.add(pos);
            }
        }
        Collections.sort(positions);

        if (!new HashSet<>(taxaNames).equals(new HashSet<>(Arrays.asList(readCount.getTaxaNames())))) {
            throw new IllegalArgumentException("Taxa names do not match read count data.");
        }

        // get column number from pos
        int[] sites = readCount.getSitesIndex();

        Map<Integer, Integer> siteToCol = new HashMap<>();
        for (int j = 0; j < sites.length; j++) {
            siteToCol.put(sites[j], j);
        }

        int[] columns = new int[positions.size()];
        for (int i = 0; i < positions.size(); i++) {
            columns[i] = siteToCol.getOrDefault(positions.get(i), -1);
        }

        // initialise a new read count matrix
        ReadCount[][] matrix = new ReadCount[taxaNames.size()][positions.size()];
        for (int i = 0; i < taxaNames.size(); i++) {
            for (int j = 0; j < columns.length; j++) {
                matrix[i][j] = readCount.getState(taxaNames.get(i), columns[j]-1); //change to 0-base index
            }
        }

        ReadCountData outputReadCount = new ReadCountData(Taxa.createTaxa(taxaNames.toArray(new String[0])),
                matrix,
                positions.stream().mapToInt(Integer::intValue).toArray());

        return new Value<>("", outputReadCount, this);
    }

    public Value<CellPosition[]> getSites() {
        return getParams().get(sitesName);
    }

    public Value<ReadCountData> getReadCount() {
        return getParams().get(readCountName);
    }
}
