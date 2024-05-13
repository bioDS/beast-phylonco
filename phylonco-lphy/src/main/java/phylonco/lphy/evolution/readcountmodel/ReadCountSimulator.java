package phylonco.lphy.evolution.readcountmodel;


import lphy.base.distribution.LogNormal;
import lphy.base.distribution.NegativeBinomial;
import lphy.base.distribution.Normal;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.coalescent.Coalescent;
import lphy.base.evolution.likelihood.PhyloCTMC;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.Value;
import phylonco.lphy.evolution.datatype.PhasedGenotype;
import phylonco.lphy.evolution.datatype.PhasedGenotypeState;
import phylonco.lphy.evolution.substitutionmodel.GT16;

public class ReadCountSimulator {
    private Value<Integer> n;
    private Value<Integer> l;
    private Value<Integer> r;
    private Value<Double> p;
    private Value<Number> M;
    private Value<Number> S;
    private Value<Double> theta;
    private Value<TimeTree> psi;
    private Value<Integer[][]> coverage;
    private Value<Alignment> data;
    private Value<Integer[][][]> readCount;
    private Multinomial multinomial = new Multinomial();
    private Value<Double> delta;
    //private Value<Number> meanDelta;
    //private Value<Number> sdDelta;
    private Value<Double> epsilon;
    //private Value<Number> meanEpsilon;
    //private Value<Number> sdEpsilon;
    private Value<Double> t;
    private Value<Number> meanT;
    private Value<Number> sdT;
    private Value<Double> v;
    private Value<Number> meanV;
    private Value<Number> sdV;
    private Value<Double> s;
    private Value<Number> meanS;
    private Value<Number> sdS;
    //private Normal normalDelta;
    //private Normal normalEpsilon;
    private LogNormal logNormalT;
    private LogNormal logNormalV;
    private LogNormal logNormalS;
    private Value<Integer[][]> alpha;
    private Multinomial multinomialAlpha;
    //Maybe generate all the Alphas one time.



    public ReadCountSimulator(Integer n, Integer l, Double M, Double S, Double delta, Double epsilon,
                              Number meanT, Number sdT, Number meanV, Number sdV, Number meanS, Number sdS){
        this.n = new Value<>("n", n);
        this.l = new Value<>("l", l);
        this.M = new Value<>("M", M);
        this.S = new Value<>("S", S);
        //this.meanDelta = new Value<>("meanDelta", meanDelta);
        //this.sdDelta = new Value<>("sdDelta", sdDelat);
        //normalDelta = new Normal(this.meanDelta, this.sdDelta);
        //delta = normalDelta.sample();
        //this.meanEpsilon = new Value<>("meanEpsilon", meanEpsilon);
        //this.sdEpsilon = new Value<>("sdEpsilon", sdEpsilon);
        //normalEpsilon = new Normal(this.meanEpsilon, this.sdEpsilon);
        //epsilon = normalEpsilon.sample();
        this.delta = new Value<>("delta", delta);
        this.epsilon = new Value<>("epsilon", epsilon);
        this.meanT = new Value<>("meanT", meanT);
        this.sdT = new Value<>("sdT", sdT);
        logNormalT = new LogNormal(this.meanT, this.sdT, null);
        t = logNormalT.sample();
        this.meanV = new Value<>("meanV", meanV);
        this.sdV = new Value<>("sdV", sdV);
        logNormalV = new LogNormal(this.meanV, this.sdV, null);
        v = logNormalV.sample();
        this.meanS = new Value<>("meanS", meanS);
        this.sdS = new Value<>("sdS", sdS);
        logNormalS = new LogNormal(this.meanS, this.sdS, null);
        Double[] proAlpha = {this.delta.value(), 1-this.delta.value()};
        Value<Double[]> proA = new Value<>("proA", proAlpha);
        Value<Integer> numberA = new Value<>("numberA", 1);
        multinomialAlpha = new Multinomial(numberA, proA);

    }



    public void simulateReadCount() {
        simulateData();
        simulateAlpha();
        simulateCoverage();
        Integer[][][] readC = new Integer[n.value()][l.value()][4];


        for (int i = 0; i < n.value(); i++) {
            for (int j = 0; j < l.value(); j++) {
                int stateIndex = data.value().getState(i, j);
                PhasedGenotypeState genotypeState = PhasedGenotype.getCanonicalState(stateIndex);
                String genotype = genotypeState.getFullName();//get genotype
                Integer cov = coverage.value()[i][j];
                Double[] proMatrix;
                Value<Integer> cover = new Value<>("cover", cov);
                multinomial.setParam("n", cover);

                if (genotype.equals("AA")){
                    Double eps = epsilon.value();
                    Double pA = 1-eps;
                    Double pC = eps/3;
                    Double pG = eps/3;
                    Double pT = eps/3;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                } else if (genotype.equals("AC") || genotype.equals("CA")){
                    Double eps = epsilon.value();
                    Double pA = 0.5-eps/3;
                    Double pC = 0.5-eps/3;
                    Double pG = eps/3;
                    Double pT = eps/3;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                } else if (genotype.equals("AG") || genotype.equals("GA")){
                    Double eps = epsilon.value();
                    Double pA = 0.5-eps/3;
                    Double pC = eps/3;
                    Double pG = 0.5-eps/3;
                    Double pT = eps/3;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                } else if (genotype.equals("AT") || genotype.equals("TA")){
                    Double eps = epsilon.value();
                    Double pA = 0.5-eps/3;
                    Double pC = eps/3;
                    Double pG = eps/3;
                    Double pT = 0.5-eps/3;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                } else if (genotype.equals("CC")){
                    Double eps = epsilon.value();
                    Double pA = eps/3;
                    Double pC = 1-eps;
                    Double pG = eps/3;
                    Double pT = eps/3;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                } else if (genotype.equals("CG") || genotype.equals("GC")){
                    Double eps = epsilon.value();
                    Double pA = eps/3;
                    Double pC = 0.5-eps/3;
                    Double pG = 0.5-eps/3;
                    Double pT = eps/3;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                } else if (genotype.equals("CT") || genotype.equals("TC")){
                    Double eps = epsilon.value();
                    Double pA = eps/3;
                    Double pC = 0.5-eps/3;
                    Double pG = eps/3;
                    Double pT = 0.5-eps/3;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                } else if (genotype.equals("GG")){
                    Double eps = epsilon.value();
                    Double pA = eps/3;
                    Double pC = eps/3;
                    Double pG = 1-eps;
                    Double pT = eps/3;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                } else if (genotype.equals("GT") || genotype.equals("TG")){
                    Double eps = epsilon.value();
                    Double pA = eps/3;
                    Double pC = eps/3;
                    Double pG = 0.5-eps/3;
                    Double pT = 0.5-eps/3;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                } else if (genotype.equals("TT")){
                    Double eps = epsilon.value();
                    Double pA = eps/3;
                    Double pC = eps/3;
                    Double pG = eps/3;
                    Double pT = 1-eps;
                    proMatrix = new Double[]{pA, pC, pG, pT};
                    Value<Double[]> proMat = new Value<>("proMatrix", proMatrix);
                    multinomial.setParam("p", proMat);
                    readC[i][j] = multinomial.sample().value();
                }

            }

        }
        readCount = new Value<>("readCount", readC);

    }


    private void simulateData() {

        Double equalRate = Double.valueOf(1.0/6);
        Double[] ratesValue = new Double[]{
                equalRate, equalRate, equalRate,
                equalRate, equalRate, equalRate};

        // Dirichlet rates for GT16 substitution model
        // rates ~ Dirichlet(conc=[1.0, 2.0, 1.0, 1.0, 2.0, 1.0]);
        Value<Double[]> gt16Rates = new Value<>("gt16Rates", ratesValue);
        Double equalFreqRate = Double.valueOf(1.0/16);
        Double[] gt16FreqsValue = new Double[]{
                equalFreqRate, equalFreqRate, equalFreqRate, equalFreqRate,
                equalFreqRate, equalFreqRate, equalFreqRate, equalFreqRate,
                equalFreqRate, equalFreqRate, equalFreqRate, equalFreqRate,
                equalFreqRate, equalFreqRate, equalFreqRate, equalFreqRate};

        // Dirichlet frequencies
        // pi ~ Dirichlet(conc=[3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0]);
        Value<Double[]> gt16Freqs = new Value<>("gt16Freqs", gt16FreqsValue);
        Value<Number> meanRate = new Value<>("meanRate", Double.valueOf(1.0));
        Value<Double[][]> Q = new GT16(gt16Rates, gt16Freqs, meanRate).apply(); // GT16

        LogNormal logNormal = new LogNormal(this.M, this.S, null);
        theta = logNormal.sample();
        Coalescent coalescent = new Coalescent(theta, this.n, null);
        psi = coalescent.sample();
        PhyloCTMC phyloCTMC = new PhyloCTMC(psi, null, null, Q,
                null, null, this.l, null, null);
        data = phyloCTMC.sample(); // alignment genotypes

    }


    private void simulateCoverage() {
        Integer[][] cov = new Integer[this.n.value()][this.l.value()];

//        Integer[] kValues = new Integer[this.n.value() * this.l.value()];
//        for (int i = 0; i < kValues.length; i++) {
//            // alpha = 1 (dropout) or 2 (no dropout)
//            // r and p in terms of mu_ij and sigma_ij
//            NegativeBinomial negativeBinomial = new NegativeBinomial(this.r, this.p);
//            kValues[i] = negativeBinomial.sample().value();
//        }
//        int index = 0;
        for (int i = 0; i < n.value(); i++) {
            this.s = logNormalS.sample();
            for (int j = 0; j < l.value(); j++) {
                Double mean = this.alpha.value()[i][j] * this.t.value() * this.s.value();
                Double variance = mean + (Math.pow((double) this.alpha.value()[i][j], 2.0) * v.value() * Math.pow((this.s.value()), 2.0));
                double p = mean / variance;
                float rFloat = (float) (Math.pow(mean, 2) / (variance - mean));
                int r = Math.round(rFloat);
                this.r = new Value<>("r", r);
                this.p = new Value<>("p", p);
                NegativeBinomial negativeBinomial = new NegativeBinomial(this.r, this.p);
                cov[i][j] = negativeBinomial.sample().value();
            }
        }

        coverage = new Value<>("coverage", cov);

    }


    private void simulateAlpha() {
        Integer[][] alp = new Integer[this.n.value()][this.l.value()];
        for (int i = 0; i < n.value(); i++) {
            for (int j = 0; j < l.value(); j++) {
                Value<Integer[]> alpha1 = multinomialAlpha.sample();
                if (alpha1.value()[0] == 1){
                    alp[i][j] = 1;
                }else alp[i][j] = 2;
            }
        }
        alpha = new Value<>("alpha", alp);
    }



    public void printData()  {
        Alignment alignment = data.value();
        int numCells = alignment.ntaxa();
        int numSites = alignment.nchar();
        System.out.println("Alignment: ");
        for (int i = 0; i < numCells; i++) {
            for (int j = 0; j < numSites; j++) {
                int stateIndex = alignment.getState(i, j);
                PhasedGenotypeState genotypeState = PhasedGenotype.getCanonicalState(stateIndex);
                String genotype = genotypeState.getFullName();
                System.out.print(genotype + "\t");
            }
            System.out.println();
        }

    }


    public void printCoverage() {
        System.out.println("Coverage: ");
        for (int i = 0; i < coverage.value().length; i++) {
            for (int j = 0; j < coverage.value()[i].length; j++) {
                System.out.print(coverage.value()[i][j] + "\t");
            }
            System.out.println();
        }
    }

    public void printAlpha() {
        System.out.println("Alpha: ");
        for (int i = 0; i < alpha.value().length; i++) {
            for (int j = 0; j < alpha.value()[i].length; j++) {
                System.out.print(alpha.value()[i][j] + "\t");
            }
            System.out.println();
        }
    }

    public Value<Integer[][]> getCoverage() {
        return this.coverage;
    }

    public Value<Alignment> getData() {
        return this.data;
    }

    public Value<Integer[][]> getAlpha() {
        return this.alpha;
    }

    public Value<Integer[][][]> getReadCount() {
        return this.readCount;
    }


}
