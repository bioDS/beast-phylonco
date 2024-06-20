package phylonco.lphy.evolution.readcountmodel;

import lphy.base.distribution.LogNormal;
import lphy.base.distribution.Multinomial;
import lphy.base.distribution.NegativeBinomial;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.Alignment;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.ParameterInfo;
import phylonco.lphy.evolution.datatype.PhasedGenotype;
import phylonco.lphy.evolution.datatype.PhasedGenotypeState;

import java.util.Map;
// lphy script
// n = 40;
// l = 200;
// t ~ LogNormal(mean, sd);
// v ~ LogNormal(mean2, sd2);
// ones ~ rep(1, n*l);
// alpha ~ Multinomial(n=1, p=[delta], replicates=n*l) + ones; // alpha_ij = 1 or 2
// sPrime ~ Lognormal(sMean, sStdDev, replicates=n);
// s ~ transpose(rep(sPrime, l));
// mean ~ alpha * rep(t, n*l) * s;
// variance ~
// p ~ // NegativeBinomial p
// r ~ // NegativeBinomial r
// cov ~ NegativeBinomial(...);
// cov ~ reshape(cov, size=[n, l]);
// data ~ PhyloCTMC();
// Reads ~ ReadCountModel(data, cov, delta, epsilon);

public class ReadCountModel implements GenerativeDistribution<ReadCountData> {

    private Value<Integer[][]> coverage;
    private Value<Alignment> data;
    //private Value<Integer[][][]> readCount;
    private Multinomial multinomial = new Multinomial();
    private Value<Double> epsilon;
//    private Value<Double> delta;
    private RandomVariable<ReadCountData> readCountData;
    private Value<Integer[][]> alpha;
//    private Multinomial multinomialAlpha;
//    private Value<Double> t;
//    private Value<Number> meanT;
//    private Value<Number> sdT;
//    private Value<Double> v;
//    private Value<Number> meanV;
//    private Value<Number> sdV;
//    private Value<Double> s;
//    private Value<Number> meanS;
//    private Value<Number> sdS;
//    private LogNormal logNormalT;
//    private LogNormal logNormalV;
//    private LogNormal logNormalS;



    public static final String dParamName = "D";
    public static final String covParamName = "coverage";
    public static final String alphaParamName = "alpha";
//    public static final String deltaParamName = "delta";
    public static final String epsilonParamName = "epsilon";
//    public static final String meanTParamName = "meanT";
//    public static final String sdTParamName = "sdT";
//    public static final String meanVParamName = "meanV";
//    public static final String sdVParamName = "sdV";
//    public static final String meanSParamName = "meanS";
//    public static final String sdSParamName = "sdS";


    public ReadCountModel(
            @ParameterInfo(name = dParamName, description = "genotype alignment.") Value<Alignment> data,
            @ParameterInfo(name = covParamName, description = "coverage.") Value<Integer[][]> coverage,
            @ParameterInfo(name = alphaParamName, description = "allelic dropout events for each cell at each site.") Value<Integer[][]> alpha,
//            @ParameterInfo(name = deltaParamName, description = "allelic dropout error probability.") Value<Double> delta,
            @ParameterInfo(name = epsilonParamName, description = "sequencing and amplification error probability.") Value<Double> epsilon
//            @ParameterInfo(name = meanTParamName, description = "sequencing and amplification error probability.") Value<Number> meanT,
//            @ParameterInfo(name = sdTParamName, description = "sequencing and amplification error probability.") Value<Number> sdT,
//            @ParameterInfo(name = meanVParamName, description = "sequencing and amplification error probability.") Value<Number> meanV,
//            @ParameterInfo(name = sdVParamName, description = "sequencing and amplification error probability.") Value<Number> sdV,
//            @ParameterInfo(name = meanSParamName, description = "sequencing and amplification error probability.") Value<Number> meanS,
//            @ParameterInfo(name = sdSParamName, description = "sequencing and amplification error probability.") Value<Number> sdS
            ) {
        super();
        this.data = data;
        this.coverage = coverage;
        this.alpha = alpha;
//        this.delta = delta;
        this.epsilon = epsilon;
//        this.meanT = meanT;
//        this.sdT = sdT;
//        this.meanV = meanV;
//        this.sdV = sdV;
//        this.meanS = meanS;
//        this.sdS = sdS;


    }

    // D ~ GT16 genotypes (Alignment)
    // t ~ ...
    // return ACGT {{1,2,3,4}, {4,5,6,7}, ... }

    @Override
    public RandomVariable<ReadCountData> sample() {
        int n = data.value().getTaxa().length();
        int l = data.value().nchar();
//        Value<Integer> r;
//        Value<Double> p;
        Double eps = epsilon.value();
        Double pA;
        Double pC;
        Double pG;
        Double pT;
        Integer[][][] readC = new Integer[n][l][4];
//        Integer[][] cov = new Integer[n][l];
//        Integer[][] alp = new Integer[n][l];
//        Double[] proAlpha = {this.delta.value(), 1-this.delta.value()};
//        Value<Double[]> proA = new Value<>("proA", proAlpha);
//        Value<Integer> numberA = new Value<>("numberA", 1);
//        multinomialAlpha = new Multinomial(numberA, proA);
//        logNormalT = new LogNormal(this.meanT, this.sdT, null);
//        logNormalV = new LogNormal(this.meanV, this.sdV, null);
//        logNormalS = new LogNormal(this.meanS, this.sdS, null);
//        t = logNormalT.sample();
//        v = logNormalV.sample();



/*        for (int i = 0; i < n; i++) {
            this.s = logNormalS.sample();
            for (int j = 0; j < l; j++) {
                Value<Integer[]> alpha1 = multinomialAlpha.sample();
                if (alpha1.value()[0] == 1){
                    alp[i][j] = 1;
                }else alp[i][j] = 2;

                Double mean = alp[i][j] * this.t.value() * this.s.value();
                Double variance = mean + (Math.pow((double) alp[i][j], 2.0) * v.value() * Math.pow((this.s.value()), 2.0));
                double pValue = mean / variance;
                float rFloat = (float) (Math.pow(mean, 2) / (variance - mean));
                int rValue = Math.round(rFloat);
                r = new Value<>("r", rValue);
                p = new Value<>("p", pValue);
                NegativeBinomial negativeBinomial = new NegativeBinomial(r, p);
                cov[i][j] = negativeBinomial.sample().value();
            }
        }
        coverage = new Value<>("coverage", cov);
        alpha = new Value<>("alpha", alp);
**/

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < l; j++) {
                int stateIndex = data.value().getState(i, j);
                PhasedGenotypeState genotypeState = PhasedGenotype.getCanonicalState(stateIndex);
                String genotype = genotypeState.getFullName();//get genotype
                Integer cov_ij = coverage.value()[i][j];
                Double[] proMatrix;
                Value<Integer> cover = new Value<>("cover", cov_ij);
                multinomial.setParam("n", cover);
                if (genotype.equals("AA")){
                    pA = 1-eps; pC = eps/3; pG = eps/3; pT = eps/3;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                } else if (genotype.equals("AC") || genotype.equals("CA")){
                    if (alpha.value()[i][j] == 1){
                        if (Math.random() < 0.5){
                            pA = 1-eps; pC = eps/3; pG = eps/3; pT = eps/3;
                            proMatrix = new Double[]{pA, pC, pG, pT};
                            Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                            multinomial.setParam("p", proMat);
                            readC[i][j] = multinomial.sample().value();
                        }else {
                            pA = eps/3; pC = 1-eps; pG = eps/3; pT = eps/3;
                            proMatrix = new Double[]{pA, pC, pG, pT};
                            Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                            multinomial.setParam("p", proMat);
                            readC[i][j] = multinomial.sample().value();
                        }
                    }else {
                        pA = 0.5 - eps / 3; pC = 0.5 - eps / 3; pG = eps / 3; pT = eps / 3;
                        proMatrix = new Double[]{pA, pC, pG, pT};
                        Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                        multinomial.setParam("p", proMat);
                        readC[i][j] = multinomial.sample().value();
                    }


                } else if (genotype.equals("AG") || genotype.equals("GA")){
                    if (alpha.value()[i][j] == 1){
                        if (Math.random() < 0.5){
                            pA = 1-eps; pC = eps/3; pG = eps/3; pT = eps/3;
                            proMatrix = new Double[]{pA, pC, pG, pT};
                            Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                            multinomial.setParam("p", proMat);
                            readC[i][j] = multinomial.sample().value();
                        }else {
                            pA = eps/3; pC = eps/3; pG = 1-eps; pT = eps/3;
                            proMatrix = new Double[]{pA, pC, pG, pT};
                            Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                            multinomial.setParam("p", proMat);
                            readC[i][j] = multinomial.sample().value();
                        }
                    }else {
                    pA = 0.5-eps/3; pC = eps/3; pG = 0.5-eps/3; pT = eps/3;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                    }


                } else if (genotype.equals("AT") || genotype.equals("TA")){
                    if (alpha.value()[i][j] == 1){
                        if (Math.random() < 0.5){
                            pA = 1-eps; pC = eps/3; pG = eps/3; pT = eps/3;
                            proMatrix = new Double[]{pA, pC, pG, pT};
                            Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                            multinomial.setParam("p", proMat);
                            readC[i][j] = multinomial.sample().value();
                        }else {
                            pA = eps/3; pC = eps/3; pG = eps/3; pT = 1-eps;
                            proMatrix = new Double[]{pA, pC, pG, pT};
                            Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                            multinomial.setParam("p", proMat);
                            readC[i][j] = multinomial.sample().value();
                        }
                    }else {
                        pA = 0.5 - eps / 3; pC = eps / 3; pG = eps / 3; pT = 0.5 - eps / 3;
                        proMatrix = new Double[]{pA, pC, pG, pT};
                        Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                        multinomial.setParam("p", proMat);
                        readC[i][j] = multinomial.sample().value();
                    }


                } else if (genotype.equals("CC")){
                    pA = eps/3; pC = 1-eps; pG = eps/3; pT = eps/3;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                }

                else if (genotype.equals("CG") || genotype.equals("GC")){
                    if (alpha.value()[i][j] == 1){
                        if (Math.random() < 0.5){
                            pA = eps/3; pC = 1-eps; pG = eps/3; pT = eps/3;
                            proMatrix = new Double[]{pA, pC, pG, pT};
                            Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                            multinomial.setParam("p", proMat);
                            readC[i][j] = multinomial.sample().value();
                        }else {
                            pA = eps/3; pC = eps/3; pG = 1-eps; pT = eps/3;
                            proMatrix = new Double[]{pA, pC, pG, pT};
                            Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                            multinomial.setParam("p", proMat);
                            readC[i][j] = multinomial.sample().value();
                        }
                    } else {
                        pA = eps / 3; pC = 0.5 - eps / 3; pG = 0.5 - eps / 3; pT = eps / 3;
                        proMatrix = new Double[]{pA, pC, pG, pT};
                        Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                        multinomial.setParam("p", proMat);
                        readC[i][j] = multinomial.sample().value();
                    }


                } else if (genotype.equals("CT") || genotype.equals("TC")){
                    if (alpha.value()[i][j] == 1){
                        if (Math.random() < 0.5){
                            pA = eps/3; pC = 1-eps; pG = eps/3; pT = eps/3;
                            proMatrix = new Double[]{pA, pC, pG, pT};
                            Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                            multinomial.setParam("p", proMat);
                            readC[i][j] = multinomial.sample().value();
                        }else {
                            pA = eps/3; pC = eps/3; pG = eps/3; pT = 1-eps;
                            proMatrix = new Double[]{pA, pC, pG, pT};
                            Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                            multinomial.setParam("p", proMat);
                            readC[i][j] = multinomial.sample().value();
                        }
                    } else {
                        pA = eps / 3; pC = 0.5 - eps / 3; pG = eps / 3; pT = 0.5 - eps / 3;
                        proMatrix = new Double[]{pA, pC, pG, pT};
                        Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                        multinomial.setParam("p", proMat);
                        readC[i][j] = multinomial.sample().value();
                    }


                } else if (genotype.equals("GG")){
                    pA = eps/3; pC = eps/3; pG = 1-eps; pT = eps/3;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                }


                else if (genotype.equals("GT") || genotype.equals("TG")){
                    if (alpha.value()[i][j] == 1){
                        if (Math.random() < 0.5){
                            pA = eps/3; pC = eps/3; pG = 1-eps; pT = eps/3;
                            proMatrix = new Double[]{pA, pC, pG, pT};
                            Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                            multinomial.setParam("p", proMat);
                            readC[i][j] = multinomial.sample().value();
                        }else {
                            pA = eps/3; pC = eps/3; pG = eps/3; pT = 1-eps;
                            proMatrix = new Double[]{pA, pC, pG, pT};
                            Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                            multinomial.setParam("p", proMat);
                            readC[i][j] = multinomial.sample().value();
                        }
                    }else {
                        pA = eps / 3; pC = eps / 3; pG = 0.5 - eps / 3; pT = 0.5 - eps / 3;
                        proMatrix = new Double[]{pA, pC, pG, pT};
                        Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                        multinomial.setParam("p", proMat);
                        readC[i][j] = multinomial.sample().value();
                    }


                } else if (genotype.equals("TT")){
                    pA = eps/3; pC = eps/3; pG = eps/3; pT = 1-eps;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                }

            }

        }
        //readCount = new Value<>("readCount", readC);
        Taxa taxa = data.value().getTaxa();
        ReadCount[][] readCountMatrix = new ReadCount[n][l];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < l; j++) {
                System.out.println("cell = " + i + ", site = " + j);
                System.out.println("coverage = " + coverage.value()[i][j]);
                int stateIndex = data.value().getState(i, j);
                PhasedGenotypeState genotypeState = PhasedGenotype.getCanonicalState(stateIndex);
                String genotype = genotypeState.getFullName();
                System.out.println("genotype = " + genotype);
                System.out.println("alpha = " + alpha.value()[i][j]);
                int[] values = new int[ReadCount.NUM_NUCLEOTIDES];
                for (int k = 0; k < ReadCount.NUM_NUCLEOTIDES; k++) {
                    values[k] = (int) readC[i][j][k];
                    System.out.print(values[k] + "\t");
                }
                System.out.println();
                readCountMatrix[i][j] = new ReadCount(values);

            }
        }
        readCountData = new RandomVariable<>("readCountData", new ReadCountData(taxa, readCountMatrix), this);

        return readCountData;
    }

    @Override
    public Map<String, Value> getParams() {
        return Map.of(
                dParamName, data,
                covParamName, coverage,
                alphaParamName, alpha,
//                deltaParamName, delta,
                epsilonParamName, epsilon
//                meanTParamName,meanT,
//                sdTParamName,sdT,
//                meanVParamName,meanV,
//                sdVParamName,sdV,
//                meanSParamName,meanS,
//                sdSParamName,sdS
        );
    }
}
