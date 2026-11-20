package phylonco.lphybeast.tobeast.values;

import beast.base.evolution.alignment.Sequence;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import phylonco.beast.evolution.alignment.ReadCountAlignment;
import phylonco.lphy.evolution.readcountmodel.ReadCountData;

/**
 * Converts LPhy read-count data into a {@link ReadCountAlignment}: one sequence per cell, sites
 * separated by ',' and the four base counts by ':'.
 *
 * <p>This is the same object BEAUti produces for an imported read-count file, so both routes into
 * BEAST generate the same shape of XML, and the read counts are the partition alignment rather than
 * a separate data object plus a genotype scaffold.</p>
 *
 * <p>The declared genotype datatype is left at the alignment's default: {@code ReadCountTreeLikelihood}
 * sets it from the GT10/GT16 substitution model, which is authoritative.</p>
 */
public class ReadCountToBEAST implements ValueToBEAST<ReadCountData, ReadCountAlignment> {

    @Override
    public ReadCountAlignment valueToBEAST(Value<ReadCountData> value, BEASTContext context) {
        ReadCountData data = value.value();
        int n = data.getTaxa().getDimension();
        int l = data.nchar();
        String[] taxaNames = data.getTaxaNames();

        ReadCountAlignment alignment = new ReadCountAlignment();
        for (int i = 0; i < n; i++) {
            StringBuilder counts = new StringBuilder();
            for (int j = 0; j < l; j++) {
                if (j > 0) {
                    counts.append(",");
                }
                counts.append(String.format("%d:%d:%d:%d",
                        data.getState(i, j).getCount("A"),
                        data.getState(i, j).getCount("C"),
                        data.getState(i, j).getCount("G"),
                        data.getState(i, j).getCount("T")));
            }
            alignment.setInputValue("sequence", new Sequence(taxaNames[i], counts.toString()));
        }

        int[] sitesIndex = data.getSitesIndex();
        StringBuilder sites = new StringBuilder();
        for (int i = 0; i < sitesIndex.length; i++) {
            if (i > 0) {
                sites.append(" ");
            }
            sites.append(sitesIndex[i] - 1); // LPhy site indices are 1-based
        }
        alignment.setInputValue("sitesIndex", sites.toString());

        alignment.initAndValidate();
        return alignment;
    }

    @Override
    public Class getValueClass() {
        return ReadCountData.class;
    }
}
