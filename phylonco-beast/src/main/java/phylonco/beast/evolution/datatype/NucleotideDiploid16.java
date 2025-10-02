package phylonco.beast.evolution.datatype;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;

@Description("Phased diploid nucleotide data type")
public class NucleotideDiploid16 extends DataType.Base implements DiploidDataType {
    int[][] x = {
            {0}, // AA - 0
            {1}, // AC - 1
            {2}, // AG - 2
            {3}, // AT - 3
            {4}, // CA - 4
            {5}, // CC - 5
            {6}, // CG - 6
            {7}, // CT - 7
            {8}, // GA - 8
            {9}, // GC - 9
            {10}, // GG - a
            {11}, // GT - b
            {12}, // TA - c
            {13}, // TC - d
            {14}, // TG - e
            {15}, // TT - f
            {1, 4}, // AC or CA - M
            {2, 8}, // AG or GA - R
            {3, 12}, // AT or TA - W
            {6, 9}, // CG or GC - S
            {7, 13}, // CT or TC - Y
            {11, 14}, // GT or TG - K
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}, // gap -
            {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15} // missing ?
    };

    public NucleotideDiploid16() {
        stateCount = 16;
        mapCodeToStateSet = x;
        codeLength = 1;
        codeMap = "0123456789ABCDEFMRWSYK" + GAP_CHAR + MISSING_CHAR;
    }

    @Override
    public String getTypeDescription() {
        return "nucleotideDiploid16";
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