package phylonco.beast.evolution.datatype;


import beast.base.core.BEASTObject;
import beast.base.core.Input;

import java.util.ArrayList;

public class ReadCount extends BEASTObject {
    int ntaxa;
    int nchar;
    int data[][][];
    String[] taxaNames;

    String newline = "\n";
    String semicolon = ",";
    String comma = ":";
    String space = " ";

    public Input<String> readCountStrInput = new Input<>("value", "A string record read counts", Input.Validate.REQUIRED);
    public Input<String> taxaNamesInput = new Input<>("taxaNames", "Taxa names of the data", Input.Validate.OPTIONAL);

    public ReadCount() {
        // do we need to fill this in?
    }

    public ReadCount(String data) {
        this.readCountStrInput = new Input("value", "A string record read counts", Input.Validate.REQUIRED);
        this.readCountStrInput.setValue(data, this);
        this.initAndValidate();
    }

    public ReadCount(int ntaxa, int nchar) {
        this.ntaxa = ntaxa;
        this.nchar = nchar;
        data = new int[ntaxa][nchar][4];
    }

    public int[] getReadCounts(int taxa, int site) {
        return data[taxa][site];
    }

    public void setReadCounts(int taxa, int site, int[] counts) {
        for (int i = 0; i < data[taxa][site].length; i++) {
            data[taxa][site][i] = counts[i];
        }
    }

    public String getTaxonName(int taxa) {return taxaNames[taxa];}


    public String getTypeDescription() {
        return "readCount";
    }

    private ArrayList<String> getTrimmed(String[] array) {
        ArrayList<String> trimmed = new ArrayList<>();
        for (String string : array) {
            String num = string.trim();
            if (num.length() > 0) {
                trimmed.add(num);
            }
        }
        return trimmed;
    }

    @Override
    public void initAndValidate() {

        // getValue from readCountStrInput
        String readCountStr = readCountStrInput.get();
        // data format A,C,G,T; A,C,G,T;
        /**
         * 1:0:0:11,0:17:0:12
         * 7:0:0:26,0:12:0:8
         */
        ArrayList<String> cellArray;
        ArrayList<String> siteArray;
        String[] cells = readCountStr.split(newline);
        cellArray = getTrimmed(cells);

        // set nchar and ntaxa
        this.ntaxa = cellArray.size();
        this.nchar = cellArray.get(0).split(semicolon).length;
        data = new int[ntaxa][nchar][4];

        int taxaIndex = 0;
        for (String cell : cellArray) {
            String[] sites = cell.split(semicolon);
            siteArray = getTrimmed(sites);
            // deal with spaces and new lines using trimmed method
            int siteIndex = 0;
            for (String site : siteArray) {
                String[] nucleotideCounts = site.split(comma);
                // convert counts to int[]
                int[] counts = new int[nucleotideCounts.length];
                for (int i = 0; i < nucleotideCounts.length; i++) {
                    counts[i] = Integer.parseInt(nucleotideCounts[i]);
                }
                // fill in counts
                setReadCounts(taxaIndex, siteIndex, counts);
                siteIndex++;
            }
            taxaIndex++;
        }

        if (taxaNamesInput.get() == null) {
            taxaNames = new String[taxaIndex];
            for (int i = 0; i < taxaIndex; i++) {
                taxaNames[i] = String.format("%d", i);
            }
        } else {
            String taxaNamesStr = taxaNamesInput.get();
            taxaNames = taxaNamesStr.split(space);
        }
    }
}



