package phylonco.lphy.evolution.readcountmodel;

import lphy.base.evolution.Mpileup;
import lphy.base.evolution.PileupSite;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.Taxon;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.lang.reflect.Array;
import java.util.*;

import static lphy.base.evolution.PileupSite.translateRead;

public class MpileupToReadCount extends DeterministicFunction<ReadCountData> {

    public MpileupToReadCount(
            @ParameterInfo(
                    name = "mpileup",
                    description = "the mpileup data converts to read count data, supports multiple chromosomes"
            )
            Value<Mpileup> mpileup) {
        if (mpileup == null) {
            throw new IllegalArgumentException("Mpileups cannot be null");
        }
        setParam("mpileup", mpileup);
    }

    @GeneratorInfo(
            name = "toReadCount",
            examples = {"mpileupToReadCount.lphy"},
            description = "Convert the mpileups to read count data (supporting multiple chromosomes)."
    )
    @Override
    public Value<ReadCountData> apply() {
        // get param
        Mpileup mpileup = getMpileups().value();

       // get all taxa names
        Set<String> taxaNamesSet = new LinkedHashSet<>();
        taxaNamesSet.addAll(mpileup.getPileupData().get(0).keySet());
        String[] taxaNames = taxaNamesSet.toArray(new String[0]);

        //generate Taxa
        Taxon[] taxonArray = new Taxon[Array.getLength(taxaNames)];
        for (int i = 0; i < taxonArray.length; i++) {
            String name = Array.get(taxaNames, i).toString();
            String species = null;
            double age = 0.0;
            taxonArray[i] = new Taxon(name, species, age);
        }

        Taxa taxa = new Taxa.Simple(taxonArray);
        ReadCount[][] readCountMatrix = new ReadCount[taxaNames.length][mpileup.getPileupData().size()];

        String[] chroms = mpileup.getChromNames();
        int[] positions = mpileup.getPositions();
        int[] refs = mpileup.getRefs();

        for (int i = 0; i < mpileup.getPileupData().size(); i++) {
            Map<String, PileupSite.CellPileupData> currentMap = mpileup.getPileupData().get(i);
            for (int j = 0; j < taxaNames.length; j++) {
                String taxon = taxaNames[j];
                PileupSite.CellPileupData cellData = currentMap.get(taxon);
                readCountMatrix[j][i] = translateReads(refs[i], cellData);
            }
        }
        ReadCountData readCountData = new ReadCountData(
                chroms, refs,
                taxa,
                readCountMatrix, positions
        );

        return new Value<>("", readCountData, this);
    }

    public Value<Mpileup> getMpileups() {
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
