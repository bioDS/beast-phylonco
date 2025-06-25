package phylonco.lphy.evolution.readcountmodel;

import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.TaxaCharacterMatrix;
import lphy.core.logger.TextFileFormatted;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;


public class ReadCountData implements TaxaCharacterMatrix<ReadCount>, TextFileFormatted {

    ReadCount[][] readCountDataMatrix; // taxa, position
    Taxa taxa;
    int[] sitesIndex;

    public ReadCountData(Taxa taxa, ReadCount[][] readCountDataMatrix) {
        this.taxa = taxa;
        this.readCountDataMatrix = readCountDataMatrix;
    }

    public ReadCountData(Taxa taxa, ReadCount[][] readCountDataMatrix, int[] sitesIndex) {
        this.taxa = taxa;
        this.readCountDataMatrix = readCountDataMatrix;
        this.sitesIndex = sitesIndex;
    }

    public void setSitesIndex(int[] sitesIndex) {
        this.sitesIndex = sitesIndex;
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
        StringBuilder sb = new StringBuilder();
        int n = getTaxa().getDimension();
        int l = nchar();
        sb.append("\n");
        //sb.append("\"");
        for (int i = 0; i < n; i++) {
            // n taxa
            sb.append(getTaxaNames()[i]);
            sb.append("\t");
            for (int j = 0; j < l; j++) {
                // n sites
                int countA = readCountDataMatrix[i][j].getCount("A");
                int countC = readCountDataMatrix[i][j].getCount("C");
                int countG = readCountDataMatrix[i][j].getCount("G");
                int countT = readCountDataMatrix[i][j].getCount("T");
                sb.append(String.format(countA + ":" + countC + ":" + countG + ":" + countT )).append("\t");
                if (j == l-1 && i != n-1){
                    sb.append("\n");
                } else if (j != l-1) {
                    //sb.append(",");
                }
            }
        }
        //sb.append("\"");
        return sb.toString();
    }

    public String[] getTaxaNames() {
        return taxa.getTaxaNames();
    }

    public int[] getSitesIndex() {
        return sitesIndex;
    }

    @Override
    public List<String> getTextForFile() {
        // TODO: not using this at the moment
        return List.of();
    }

    @Override
    public String getFileType() {
        return "_readcount.log";
    }

    public void writeToFile(BufferedWriter writer) {
        try {
            int n = getTaxa().getDimension();
            int l = nchar();
            for (int i = 0; i < n; i++) {
                // n taxa
                for (int j = 0; j < l; j++) {
                    // n sites
                    int countA = readCountDataMatrix[i][j].getCount("A");
                    int countC = readCountDataMatrix[i][j].getCount("C");
                    int countG = readCountDataMatrix[i][j].getCount("G");
                    int countT = readCountDataMatrix[i][j].getCount("T");
                    String readCountForSiteForCell = countA + ":" + countC + ":" + countG + ":" + countT;
                    if (j == l-1 && i != n-1){
                        readCountForSiteForCell += "\n";
                    } else if (j != l-1) {
                        readCountForSiteForCell += ",";
                    }
                    // write out read count for site and cell
                    writer.write(readCountForSiteForCell);
                }
            }

            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
