package phylonco.lphy.evolution.copynumbermodel;

import org.junit.Test;
import phylonco.lphy.evolution.alignment.HomozygousAlignmentDistribution;

public class NestedBDTest {


    /**
     *
     */
    @Test
    public void simulateCopiesOnBranchBinTest() {
        NestedBD model = new NestedBD();
        // work out what the expected value is (mean for number of copies)
        int result = model.simulateCopiesOnBranchBin(1.0, 1.5, 2);

    }

}
