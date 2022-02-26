package phylonco.beast.evolution.operators;

import beast.core.parameter.RealParameter;
import beast.evolution.substitutionmodel.Frequencies;
import org.junit.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.Assert.assertArrayEquals;

public class BlockPairSwapOperatorTest {

    private static double DELTA = 1e-10;

    @Test
    public void testPairSwap() {
        Double[] values = {0.1, 0.2, 0.3, 0.05, 0.15, 0.2};
        Frequencies freq = new Frequencies();
        RealParameter f = new RealParameter(values);
        freq.initByName("frequencies", f);
        f.initAndValidate(); // extra init checking should not break code

        String pairs = "0 1 1 2 2 0";

        BlockPairSwapOperator operator = new BlockPairSwapOperator();
        operator.initByName(
                "parameter", f,
                "weight", 1.0,
                "pairs", pairs);
        operator.initAndValidate(); // extra init checking should not break code
        operator.proposal();

        // expected frequencies after swap
        double[] expected = {0.1, 0.3, 0.2, 0.15, 0.05, 0.2};
        double[] observed = Arrays.stream(f.getValues()).mapToDouble(Double::doubleValue).toArray();

        assertArrayEquals(expected, observed, DELTA);
    }

    @Test
    public void testPairSwap2() {
        Double[] values = {0.1, 0.2, 0.3, 0.05, 0.15, 0.2};
        Frequencies freq = new Frequencies();
        RealParameter f = new RealParameter(values);
        freq.initByName("frequencies", f);

        String pairs = "2 1 0 2 1 0";

        BlockPairSwapOperator operator = new BlockPairSwapOperator();
        operator.initByName(
                "parameter", f,
                "weight", 1.0,
                "pairs", pairs);
        operator.proposal();

        // expected frequencies after swap
        double[] expected = {0.05, 0.15, 0.3, 0.1, 0.2, 0.2};
        double[] observed = Arrays.stream(f.getValues()).mapToDouble(Double::doubleValue).toArray();

        assertArrayEquals(expected, observed, DELTA);

    }

    @Test
    public void testGT16Pairs() {
        // convert swap pair vectors for GT16 unphased data
        // swap AC and CA
        String[] keys = {"AC", "AG", "AT", "GC", "TC"};
        convertPairFormat(keys);

        // swap AG and GA
        keys = new String[] {"AG", "AC", "AT", "GC", "GT"};
        convertPairFormat(keys);

        // swap AT and TA
        keys = new String[] {"AT", "AC", "AG", "CT", "GT"};
        convertPairFormat(keys);

        // swap CG and GC
        keys = new String[] {"CG", "CA", "CT", "AG", "TG"};
        convertPairFormat(keys);

        // swap CT and TC

        // swap GT and TG
    }

    private void convertPairFormat(String[] keys) {
        List<String> names = Arrays.asList(new String[]{
                "AA", "AC", "AG", "AT",
                "CA", "CC", "CG", "CT",
                "GA", "GC", "GG", "GT",
                "TA", "TC", "TG", "TT"
        });

        int[] values = new int[names.size()];

        int count = 1;

        for (String key: keys) {
            int keyIndex = names.indexOf(key);
            StringBuilder str = new StringBuilder();
            str.append(key);
            String reversed = str.reverse().toString();
            int reverseIndex = names.indexOf(reversed);
            values[keyIndex] = count;
            values[reverseIndex] = count;
            count++;
        }

        System.out.println("keys: " + Arrays.asList(keys).toString());

        for (String s: names) {
            System.out.print(s + "\t");
        }
        System.out.println();
        for (int i: values) {
            System.out.print(i + "\t");
        }
        System.out.println();
    }
}
