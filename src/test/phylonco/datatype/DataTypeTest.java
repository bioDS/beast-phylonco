package phylonco.datatype;

import beast.evolution.datatype.DataType;
import beast.evolution.datatype.NucleotideDiploid10;
import beast.evolution.datatype.NucleotideDiploid16;
import beast.evolution.datatype.Ternary;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import static org.junit.Assert.assertEquals;

@RunWith(Parameterized.class)
public class DataTypeTest {

        private DataType d;
        private String s;
        private List<Integer> c;

        public DataTypeTest(DataType dataType, String sequence, List<Integer> codes) {
            d = dataType;
            s = sequence;
            c = codes;
        }

        @Parameterized.Parameters
        public static Collection<Object[]> typestrings() {
            Ternary dataTernary = new Ternary();
            NucleotideDiploid10 dataNucleotide10 = new NucleotideDiploid10();
            NucleotideDiploid16 dataNucleotide16 = new NucleotideDiploid16();

            return Arrays.asList(new Object[][] {
                    { dataTernary, "012-?102", Arrays.asList(new Integer[] { 0, 1, 2, 3, 4, 1, 0, 2 }) },

                    { dataNucleotide10, "AMRWCSYGKT-?", Arrays.asList(
                            new Integer[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11 }) },

                    { dataNucleotide16, "0123456789ABCDEFMRWSYK-?", Arrays.asList(
                            new Integer[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23 }) }
            });
        }

        @Test
        public void testStingToEncoding() {
            assertEquals(c, d.stringToEncoding(s));
        }

        @Test
        public void testRoundTrip() {
            assertEquals(s, d.encodingToString(d.stringToEncoding(s)));
        }
}

