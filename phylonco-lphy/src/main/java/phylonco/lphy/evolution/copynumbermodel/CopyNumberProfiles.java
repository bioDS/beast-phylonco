package phylonco.lphy.evolution.copynumbermodel;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.TaxaCharacterMatrix;

import java.util.Arrays;

public class CopyNumberProfiles implements TaxaCharacterMatrix<Integer> {

    private int ntaxa;
    private int nBins;

    private Taxa taxa;

    private int[][] data;

    public CopyNumberProfiles(Taxa t, int nBins) {
        this.taxa = t;

        // dimensions
        this.nBins = nBins;
        this.ntaxa = t.ntaxa();

        this.data = new int[ntaxa][nBins];
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
}
