package test.beast.evolution.datatype;

import beast.evolution.datatype.*;
import beast.evolution.datatype.DataType.Base;

import org.junit.Test;
import static org.junit.Assert.assertEquals;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

@RunWith(Parameterized.class)
public class DataTypeDeEncodeTest {
    private Base dataType;
    private String sequence;
    private List<Integer> codes;

    public DataTypeDeEncodeTest(Base dataType, String sequence, List<Integer> codes) {
        this.dataType = dataType;
        this.sequence = sequence;
        this.codes = codes;
    }


    @Parameterized.Parameters
    public static Collection<Object[]> typestrings() {
        Base NucMeth = new NucleotideMethylation();
        Base NucDip = new NucleotideDiploid();
        Base QuinaryType = new Quinary();
        Base TernaryType = new Ternary();

        return Arrays.asList(new Object[][] {
                { NucMeth, "ACGTPJO1WNX-?", Arrays.asList( new Integer[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}) },
                { NucDip, "AMRWCSYGKT?", Arrays.asList( new Integer[] {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10})},
                { QuinaryType, "012345?", Arrays.asList( new Integer[] {0, 1, 2, 3, 4, 5, 6})},
                { TernaryType, "012?", Arrays.asList( new Integer[] {0, 1, 2, 3})}
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