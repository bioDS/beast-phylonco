package phylonco.beast.evolution.datatype;


public class ReadCount {
    int ntaxa;
    int nchar;
    int data[][][];

    public ReadCount(int ntaxa, int nchar) {
        this.ntaxa = ntaxa;
        this.nchar = nchar;
        data = new int[ntaxa][nchar][4];
    }

    public int[] getReadCounts(int taxa, int site) {
        return data[taxa][site];
    }

    public void setReadCounts(int taxa, int site, int[] counts) {
        data[taxa][site] = counts;
    }


    public String getTypeDescription() {
        return "readCount";
    }
}



