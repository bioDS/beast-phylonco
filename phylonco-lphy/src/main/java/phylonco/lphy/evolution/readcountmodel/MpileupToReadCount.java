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
                    description = "the mpileup converts to read count data, supports multiple chromosomes"
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
            description = "Convert the mpileups to read count data (supporting multiple chromosomes)."
    )
    @Override
    public Value<ReadCountData> apply() {
        // get param
        List<Mpileup> mpileups = getMpileups().value();

       // group mpileups by chom name
        Map<String, List<Mpileup>> chromToMpileups = new LinkedHashMap<>();
        for (Mpileup mp : mpileups) {
            chromToMpileups.computeIfAbsent(mp.getChromName(), k -> new ArrayList<>()).add(mp);
        }

       // get all taxa names
        Set<String> taxaNamesSet = new LinkedHashSet<>();
        for (Mpileup mp : mpileups) {
            taxaNamesSet.addAll(mp.getPileupData().keySet());
        }

        List<String> taxaNames = new ArrayList<>(taxaNamesSet);
        Collections.sort(taxaNames);

        List<Integer> allPositions = new ArrayList<>();
        List<String> allChroms = new ArrayList<>();
        List<Integer> allRefs = new ArrayList<>();
        Map<Integer, Map<String, ReadCount>> posToCellReads = new LinkedHashMap<>();

        for (String chrom : chromToMpileups.keySet()) {
            List<Mpileup> chromMpileups = chromToMpileups.get(chrom);

            for (Mpileup mp : chromMpileups) {
                int pos = mp.getPosition();
                int ref = mp.getRef();
                allPositions.add(pos);
                allChroms.add(chrom);
                allRefs.add(ref);

                Map<String, ReadCount> cellMap = new LinkedHashMap<>();
                for (Map.Entry<String, PileupSite.CellPileupData> entry : mp.getPileupData().entrySet()) {
                    String cell = entry.getKey();
                    PileupSite.CellPileupData data = entry.getValue();
                    ReadCount rc = translateReads(ref, data);
                    cellMap.put(cell, rc);
                }
                posToCellReads.put(pos, cellMap);
            }
        }
        // sort it by chrom and then pos
        List<Integer> sortedIndices = new ArrayList<>();
        for (int i = 0; i < allPositions.size(); i++) sortedIndices.add(i);
        sortedIndices.sort(Comparator
                .comparing((Integer i) -> allChroms.get(i))
                .thenComparing(i -> allPositions.get(i)));

        List<Integer> positions = new ArrayList<>();
        List<String> chromArrayList = new ArrayList<>();
        List<Integer> refArrayList = new ArrayList<>();
        for (int i : sortedIndices) {
            positions.add(allPositions.get(i));
            chromArrayList.add(allChroms.get(i));
            refArrayList.add(allRefs.get(i));
        }

        int nSites = positions.size();
        int nTaxa = taxaNames.size();

        ReadCount[][] matrix = new ReadCount[nTaxa][nSites];
        for (int j = 0; j < nTaxa; j++) {
            String cell = taxaNames.get(j);
            for (int i = 0; i < nSites; i++) {
                int pos = positions.get(i);
                Map<String, ReadCount> siteMap = posToCellReads.get(pos);
                ReadCount rc = (siteMap != null) ? siteMap.get(cell) : null;
                if (rc == null) rc = new ReadCount(new int[]{0, 0, 0, 0});
                matrix[j][i] = rc;
            }
        }

        int[] siteIndex = new int[nSites];
        for (int i = 0; i < nSites; i++) {
            siteIndex[i] = positions.get(i) - 1; // change to 0-based
        }
        String[] chromArray = chromArrayList.toArray(new String[0]);
        int[] refArray = refArrayList.stream().mapToInt(Integer::intValue).toArray();

        ReadCountData readCountData = new ReadCountData(
                chromArray, refArray,
                Taxa.createTaxa(taxaNames.toArray(new String[0])),
                matrix, siteIndex
        );

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
