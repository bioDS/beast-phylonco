package phylonco.lphy.evolution.readcountmodel;

import lphy.base.evolution.Mpileup;
import lphy.base.evolution.PileupSite;
import lphy.core.model.Value;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static org.junit.Assert.assertEquals;

public class MpileupToRcTest {
    @Test
    void toRCtest() {
        List<Map<String, PileupSite.CellPileupData>> pileupData = new ArrayList<>();
        String[] chroms = new String[]{"chr1", "chr1"};
        int[] pos = new int[]{5, 10};
        int[] refs = new int[]{2, 0};
        PileupSite.CellPileupData data1 = new PileupSite.CellPileupData(4, "..TA", "~~~~");
        Map<String, PileupSite.CellPileupData> map1 = new HashMap<>();
        map1.put("cell1", data1);
        //Mpileup pileup1 = new Mpileup("chr1", 5, 2, map1);
        pileupData.add(map1);

        PileupSite.CellPileupData data2 = new PileupSite.CellPileupData(6, "..GC..", "~~~~");
        Map<String, PileupSite.CellPileupData> map2 = new HashMap<>();
        map2.put("cell1", data2);
        //Mpileup pileup2 = new Mpileup("chr1", 10, 0, map2);
        pileupData.add(map2);

        Mpileup m = new Mpileup(chroms, pos, refs, pileupData);

        MpileupToReadCount test = new MpileupToReadCount(new Value<>("", m));
        ReadCountData rc = test.apply().value();

        assertEquals((Integer) 2, rc.nchar());
        assertEquals(5, rc.getSitesIndex()[0]);
        assertEquals(10, rc.getSitesIndex()[1]);
        assertEquals(2, rc.getRefIndex()[0]);
        assertEquals(0, rc.getRefIndex()[1]);
        for (String name: rc.getChromNames()){
            assertEquals("chr1", name);
        }

        ReadCount[][] matrix = rc.getReadCountDataMatrix();
        int[] test1 = matrix[0][0].getReadCounts();
        assertEquals(1, test1[0]);
        assertEquals(0, test1[1]);
        assertEquals(2, test1[2]);
        assertEquals(1, test1[3]);

        int[] test2 = matrix[0][1].getReadCounts();
        assertEquals(4, test2[0]);
        assertEquals(1, test2[1]);
        assertEquals(1, test2[2]);
        assertEquals(0, test2[3]);
    }
}
