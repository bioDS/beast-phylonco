package phylonco.lphy.evolution.readcountmodel;

import lphy.base.evolution.Mpileup;
import lphy.base.evolution.PileupSite;
import lphy.base.evolution.Taxa;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.*;

import static lphy.base.evolution.PileupSite.translateRead;

public class MpileupToReadCount extends DeterministicFunction<ReadCountData> {

    public MpileupToReadCount(
            @ParameterInfo(
                    name = "mpileup",
                    description = "the mpileup converts to read count data, only support one chromosome in one mpileup"
            )
            Value<List<Mpileup>> mpileups) {
        if (mpileups == null) {
            throw new IllegalArgumentException("Mpileups cannot be null");
        }
        setParam("mpileup", mpileups);
    }

    @GeneratorInfo(
            name = "toReadCount",
            examples = {"mpileupToReadCount.lphy"},
            description = "Convert the mpileups to read count data (single chromosome only)."
    )
    @Override
    public Value<ReadCountData> apply() {
        List<Mpileup> mpileups = getMpileups().value();

        // ---------- check chromosome ----------
        Set<String> chromNames = new HashSet<>();
        for (Mpileup m : mpileups) {
            chromNames.add(m.getChromName());
        }
        if (chromNames.size() != 1) {
            throw new IllegalArgumentException("Mpileups must contain exactly one chromosome, found: " + chromNames);
        }


        Set<String> taxaNamesSet = new LinkedHashSet<>();
        List<Integer> positions = new ArrayList<>();
        Map<Integer, Map<String, ReadCount>> posToCellReads = new LinkedHashMap<>();

        // get taxa names and positions
        for (Mpileup mp : mpileups) {
            int pos = mp.getPosition();
            int ref = mp.getRef();
            positions.add(pos);

            Map<String, ReadCount> cellMap = new LinkedHashMap<>();

            for (Map.Entry<String, PileupSite.CellPileupData> entry : mp.getPileupData().entrySet()) {
                String cell = entry.getKey();
                PileupSite.CellPileupData data = entry.getValue();
                taxaNamesSet.add(cell);

                ReadCount rc = translateReads(ref, data);
                cellMap.put(cell, rc);
            }
            posToCellReads.put(pos, cellMap);
        }

        List<String> taxaNames = new ArrayList<>(taxaNamesSet);
        Collections.sort(positions);

        int nSites = positions.size();
        int nTaxa = taxaNames.size();

        // build rc matrix
        ReadCount[][] matrix = new ReadCount[nTaxa][nSites];
        for (int j = 0; j < nTaxa; j++) {
            String cell = taxaNames.get(j);
            for (int i = 0; i < nSites; i++) {
                int pos = positions.get(i);
                Map<String, ReadCount> siteMap = posToCellReads.get(pos);

                ReadCount rc = null;
                if (siteMap != null) {
                    rc = siteMap.get(cell);
                }
                if (rc == null) {
                    rc = new ReadCount(new int[]{0, 0, 0, 0});
                }
                matrix[j][i] = rc;
            }
        }

        ReadCountData readCountData = new ReadCountData(Taxa.createTaxa(taxaNames.toArray(new String[0])),
                matrix,
                positions.stream().mapToInt(Integer::intValue).toArray());

        return new Value<>("", readCountData, this);
    }

    public Value<List<Mpileup>> getMpileups() {
        return getParams().get("mpileup");
    }

    public static ReadCount translateReads(int ref, PileupSite.CellPileupData data) {
        String read = translateRead(ref, data.reads());
        String[] parts = read.split(":");
        int[] counts = new int[4];
        for (int i = 0; i < counts.length; i++) {
            counts[i] = Integer.parseInt(parts[i].substring(1));
        }
        return new ReadCount(counts);
    }
}
