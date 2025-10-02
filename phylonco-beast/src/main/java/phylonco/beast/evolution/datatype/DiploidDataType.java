package phylonco.beast.evolution.datatype;

public interface DiploidDataType {

    int getIndex(String genotype);

    int[] getIndices(String[] genotypes);

}
