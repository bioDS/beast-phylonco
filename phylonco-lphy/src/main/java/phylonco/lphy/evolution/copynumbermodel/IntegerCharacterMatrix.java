package phylonco.lphy.evolution.copynumbermodel;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.TaxaCharacterMatrix;
import lphy.core.model.annotation.GeneratorCategory;
import lphy.core.model.annotation.MethodInfo;

/**
 * A character matrix for integer-valued traits across taxa.
 *
 * <p>Stores discrete integer states for multiple taxa across multiple sites or bins.</p>
 *
 * <p>This is analogous to {@link lphy.base.evolution.alignment.Alignment} for
 * molecular sequences, but designed for integer-valued data.</p>
 */
public class IntegerCharacterMatrix implements TaxaCharacterMatrix<Integer> {

    private int ntaxa;
    private int nBins;
    private Taxa taxa;
    private int[][] data;


    public IntegerCharacterMatrix(Taxa t, int nBins) {
        this.taxa = t;

        // dimensions
        this.nBins = nBins;
        this.ntaxa = t.ntaxa();
        this.data = new int[ntaxa][nBins];
    }

    @MethodInfo(description = "the taxa of the copy number matrix.", narrativeName = "list of taxa",
            category = GeneratorCategory.TAXA_ALIGNMENT)
    public Taxa taxa() {
        return getTaxa();
    }

    @Override
    public Integer getState(String taxon, int characterColumn) {
        String[] taxaNames = getTaxa().getTaxaNames();
        int taxaIndex = -1;
        // Arrays.stream(taxaNames).toList().indexOf(taxon);
        for (int i = 0; i < taxaNames.length; i++) {
            if (taxaNames[i].equals(taxon)) {
                // this is the taxa index
                taxaIndex = i;
            }
        }
        return this.data[taxaIndex][characterColumn];
    }

    @Override
    public void setState(int taxon, int position, Integer state) {
        data[taxon][position] = state;
    }

    @Override
    public Class getComponentType() {
        return Integer.class;
    }

    @Override
    public String toJSON() {
        return "";
    }

    @Override
    public Taxa getTaxa() {
        return this.taxa;
    }

    @Override
    public Integer nchar() {
        return this.nBins;
    }

    @Override
    public int getDimension() {
        return 2;
    }

    @Override
    public String toString() {
        StringBuilder res = new StringBuilder();
        for (int i = 0; i < ntaxa; i++) {
            for (int j = 0; j < nBins; j++) {
                res.append(data[i][j]);
                // Add separator between values
                if (j < nBins - 1) res.append(" ");
            }
            // Add newline between taxa
            res.append("\n");
        }
        return res.toString();
    }
}
