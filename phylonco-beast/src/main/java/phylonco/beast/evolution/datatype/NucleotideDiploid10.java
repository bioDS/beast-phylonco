package phylonco.beast.evolution.datatype;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;

@Description("Unphased diploid nucleotide data type")
public class NucleotideDiploid10 extends DataType.Base implements DiploidDataType{
    int[][] x = {
            {0}, // AA - A
            {1}, // AC - M
            {2}, // AG - R
            {3}, // AT - W
            {4}, // CC - C
            {5}, // CG - S
            {6}, // CT - Y
            {7}, // GG - G
            {8}, // GT - K
            {9}, // TT - T
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, // gap -
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9} // missing ?
    };

    public NucleotideDiploid10() {
        stateCount = 10;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "AMRWCSYGKT" + GAP_CHAR + MISSING_CHAR;
    }

    @Override
    public String getTypeDescription() {
        return "nucleotideDiploid10";
    }

    @Override
    public int getIndex(String genotype) {
        switch (genotype) {
            case "AA":
                return 0;
            case "A_":
                return getIndex("AA");
            case "_A":
                return getIndex("AA");
            case "AC":
                return 1;
            case "AG":
                return 2;
            case "AT":
                return 3;
            case "CC":
                return 4;
            case "C_":
                return getIndex("CC");
            case "_C":
                return getIndex("CC");
            case "CG":
                return 5;
            case "CT":
                return 6;
            case "GG":
                return 7;
            case "G_":
                return getIndex("GG");
            case "_G":
                return getIndex("GG");
            case "GT":
                return 8;
            case "TT":
                return 9;
            case "T_":
                return getIndex("TT");
            case "_T":
                return getIndex("TT");
        }
        throw new IllegalArgumentException("Unknown genotype: " + genotype);
    }

    @Override
    public int[] getIndices(String[] genotypes) {
        int[] indices = new int[genotypes.length];
        for (int i = 0; i < genotypes.length; i++) {
            indices[i] = getIndex(genotypes[i]);
        }
        return indices;
    }
}
