package beast.evolution.datatype;

import beast.base.evolution.datatype.DataType.Base;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import static org.junit.Assert.assertEquals;

@RunWith(Parameterized.class)
public class DataTypeTest {
    private Base dataType;
    private String sequence;
    private List<Integer> codes;

    public DataTypeTest(Base dataType, String sequence, List<Integer> codes) {
        this.dataType = dataType;
        this.sequence = sequence;
        this.codes = codes;
    }


    @Parameterized.Parameters
    public static Collection<Object[]> typestrings() {
        NucleotideMethylation NucMeth = new NucleotideMethylation();
        NucleotideDiploid10 dataNucleotide10 = new NucleotideDiploid10();
        NucleotideDiploid16 dataNucleotide16 = new NucleotideDiploid16();
        Ternary TernaryType = new Ternary();

        return Arrays.asList(new Object[][] {
                { NucMeth, "ACGTPJO1WNX-?", Arrays.asList(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) },
                { dataNucleotide10, "AMRWCSYGKT-?", Arrays.asList(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11) },
                { dataNucleotide16, "0123456789ABCDEFMRWSYK-?", Arrays.asList(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23) },
                { TernaryType, "012-?", Arrays.asList(0, 1, 2, 3, 4) }
        });
    }


    @Test
    public void testStringToEncoding() {
        assertEquals(codes, dataType.stringToEncoding(sequence));
    }


    @Test
    public void testRoundTrip() {
        assertEquals(sequence, dataType.encodingToString(dataType.stringToEncoding(sequence)));
    }
}