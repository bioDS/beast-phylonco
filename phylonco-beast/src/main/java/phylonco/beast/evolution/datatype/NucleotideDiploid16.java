package phylonco.beast.evolution.datatype;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;

import java.util.Map;

import static java.util.Map.entry;

@Description("Phased diploid nucleotide data type")
public class NucleotideDiploid16 extends DataType.Base implements DiploidDataType {

    /**
     * indices need to match data in int[][] x
     */
    private static final Map<Integer, String> integerToGenotypeMap = Map.ofEntries(
            entry(0, "AA"),
            entry(1, "AC"),
            entry(2, "AG"),
            entry(3, "AT"),
            entry(4, "CA"),
            entry(5, "CC"),
            entry(6, "CG"),
            entry(7, "CT"),
            entry(8, "GA"),
            entry(9, "GC"),
            entry(10, "GG"),
            entry(11, "GT"),
            entry(12, "TA"),
            entry(13, "TC"),
            entry(14, "TG"),
            entry(15, "TT")
    );

    /**
     * data indices in int[][] x need to match getIndex(String genotype) and integerToGenotypeMap
     */
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

    /**
     * indices need to match data indices in int[][] x and integerToGenotypeMap
     */
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
            case "CA":
                return 4;
            case "CC":
                return 5;
            case "C_":
                return getIndex("CC");
            case "_C":
                return getIndex("CC");
            case "CG":
                return 6;
            case "CT":
                return 7;
            case "GA":
                return 8;
            case "GC":
                return 9;
            case "GG":
                return 10;
            case "G_":
                return getIndex("GG");
            case "_G":
                return getIndex("GG");
            case "GT":
                return 11;
            case "TA":
                return 12;
            case "TC":
                return 13;
            case "TG":
                return 14;
            case "TT":
                return 15;
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

    @Override
    public String getGenotype(int index) {
        return integerToGenotypeMap.get(index);
    }
}