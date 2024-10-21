package phylonco.lphy.evolution.readcountmodel;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.TaxaCharacterMatrix;



public class ReadCountData implements TaxaCharacterMatrix<ReadCount> {

    ReadCount[][] readCountDataMatrix; // taxa, position
    Taxa taxa;

    public ReadCountData(Taxa taxa, ReadCount[][] readCountDataMatrix) {
        this.taxa = taxa;
        this.readCountDataMatrix = readCountDataMatrix;
    }

    @Override
    public ReadCount getState(String taxon, int column) {
        return readCountDataMatrix[taxa.indexOfTaxon(taxon)][column];
    }

    @Override
    public void setState(int taxon, int position, ReadCount state) {
        // Change it if setState requires
        throw new UnsupportedOperationException("Not supported in ReadCountData yet !");
    }

    public ReadCount getState(int taxonIndex, int column) {
        return readCountDataMatrix[taxonIndex][column];
    }

    @Override
    public Class getComponentType() {
        return ReadCount.class;
    }

    @Override
    public ReadCount[] getCharacterSequence(String taxon) {
        return TaxaCharacterMatrix.super.getCharacterSequence(taxon);
    }

    @Override
    public String toJSON() {
        return "";
    }

    @Override
    public Taxa getTaxa() {
        return taxa;
    }

    @Override
    public Integer nchar() {
        return readCountDataMatrix[0].length;
    }

    @Override
    public int getDimension() {
        return nchar() * taxa.getDimension() * ReadCount.NUM_NUCLEOTIDES;
    }

    @Override
    public String toString() {
        String string = "\n";
        int n = getTaxa().getDimension();;
        int l = nchar();
        for (int i = 0; i < n; i++) {
            // n taxa
            for (int j = 0; j < l; j++) {
                // n sites
                int countA = readCountDataMatrix[i][j].getCount("A");
                int countC = readCountDataMatrix[i][j].getCount("C");
                int countG = readCountDataMatrix[i][j].getCount("G");
                int countT = readCountDataMatrix[i][j].getCount("T");
                string += String.format("A: %d, C: %d, G: %d, T: %d; \t", countA, countC, countG, countT);
            }
            string += "\n";
        }
        return string;
    }
}
