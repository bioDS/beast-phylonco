package phylonco.lphy.evolution.readcountmodel;

import lphy.base.evolution.alignment.Alignment;
import lphy.core.model.Value;
import phylonco.lphy.evolution.datatype.PhasedGenotype;
import phylonco.lphy.evolution.datatype.PhasedGenotypeState;

public class testReadCount {
    public static void main(String[] args) {

//        double mean = 5.0;
//        double variance = 5.0;
//        double p = mean / (variance * variance);
//        float rFloat = (float) (Math.pow(mean, 2) / (Math.pow(variance, 2) - mean));
//        int r = Math.round(rFloat);

        ReadCountSimulator readCountSimulator = new ReadCountSimulator(3,20,1.0,1.0, 0.163, 0.002,20, 10,2,1,2.5,1);
        readCountSimulator.simulateReadCount();
        Value<Integer[][][]> readCount = readCountSimulator.getReadCount();
        //readCountSimulator.printData();
        //readCountSimulator.printCoverage();
        Value<Integer[][]> coverage = readCountSimulator.getCoverage();
        Value<Integer[][]> alpha = readCountSimulator.getAlpha();
        Value<Alignment> data = readCountSimulator.getData();

        readCountSimulator.printData();
        readCountSimulator.printCoverage();
        readCountSimulator.printAlpha();

        for (int i = 1; i < coverage.value().length; i++) {
            for (int j = 0; j < coverage.value()[i].length; j++) {
                int stateIndex = data.value().getState(i, j);
                PhasedGenotypeState genotypeState = PhasedGenotype.getCanonicalState(stateIndex);
                String genotype = genotypeState.getFullName();
                System.out.println("Cell: " + (i+1) + ", Site: " + (j+1));
                System.out.println("Genotype: " + genotype + ", Alpha: " + alpha.value()[i][j] + ", Coverage: " + coverage.value()[i][j]);
                System.out.println("A: " + readCount.value()[i][j][0] + ", C: " + readCount.value()[i][j][1] + ", G: " + readCount.value()[i][j][2] + ", T: " + readCount.value()[i][j][3]);
                System.out.println();
            }
        }



    }
}
