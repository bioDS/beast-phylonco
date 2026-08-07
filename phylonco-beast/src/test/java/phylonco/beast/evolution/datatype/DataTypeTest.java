package phylonco.beast.evolution.datatype;

import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.DataType.Base;
import beast.pkgmgmt.BEASTClassLoader;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.List;
import java.util.Set;
import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class DataTypeTest {

    static Stream<Arguments> typestrings() {
        NucleotideMethylation NucMeth = new NucleotideMethylation();
        NucleotideDiploid10 dataNucleotide10 = new NucleotideDiploid10();
        NucleotideDiploid16 dataNucleotide16 = new NucleotideDiploid16();
        Ternary TernaryType = new Ternary();

        return Stream.of(
                Arguments.of(NucMeth, "ACGTPJO1WNX-?", List.of(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)),
                Arguments.of(dataNucleotide10, "AMRWCSYGKT-?", List.of(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)),
                Arguments.of(dataNucleotide16, "0123456789ABCDEFMRWSYK-?", List.of(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)),
                Arguments.of(TernaryType, "012-?", List.of(0, 1, 2, 3, 4))
        );
    }

    @Test
    public void testDataTypes() {
        Set<String> dataTypes = BEASTClassLoader.loadService(DataType.class);
        assertTrue(dataTypes.contains("phylonco.beast.evolution.datatype.NucleotideDiploid10"), "NucleotideDiploid10");
        assertTrue(dataTypes.contains("phylonco.beast.evolution.datatype.NucleotideDiploid16"), "NucleotideDiploid16");
    }

    @ParameterizedTest
    @MethodSource("typestrings")
    public void testStringToEncoding(Base dataType, String sequence, List<Integer> codes) {
        assertEquals(codes, dataType.stringToEncoding(sequence));
    }

    @ParameterizedTest
    @MethodSource("typestrings")
    public void testRoundTrip(Base dataType, String sequence, List<Integer> codes) {
        assertEquals(sequence, dataType.encodingToString(dataType.stringToEncoding(sequence)));
    }
}
