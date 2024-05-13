package phylonco.lphy.evolution.readcountmodel;

public class ReadCount {

    public final static int NUM_NUCLEOTIDES = 4;

    private final String nucleotides = "ACGT";

    private int[] readCounts;

    public ReadCount(int[] readCounts) {
        this.readCounts = readCounts;
    }

    public int[] getReadCounts() {
        return readCounts;
    }

    public int getCount(String nucleotide) {
        String capitalizedNucleotide = nucleotide.toUpperCase();
        int index = nucleotides.indexOf(capitalizedNucleotide);
        if (index != -1) {
            return readCounts[index];
        } else {
            return 0;
        }
    }

    public int getNumStates() {
        return NUM_NUCLEOTIDES;
    }


}
