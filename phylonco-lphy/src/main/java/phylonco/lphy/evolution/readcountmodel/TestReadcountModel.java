package phylonco.lphy.evolution.readcountmodel;

import lphy.base.distribution.LogNormal;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.coalescent.Coalescent;
import lphy.base.evolution.likelihood.PhyloCTMC;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.Value;
import phylonco.lphy.evolution.substitutionmodel.GT16;

public class TestReadcountModel {
    public static void main(String[] args) {
        Value<Integer> n = new Value<>("n", 5);
        Value<Integer> l = new Value<>("l", 10);
        Value<Double> theta;
        Value<TimeTree> psi;
        Value<Alignment> data;
        Value<Number> M = new Value<>("M",1.0);
        Value<Number> S = new Value<>("S",1.0);
        Value<Double> delta = new Value<>("delta", 0.24);
        Value<Double> epsilon = new Value<>("epsilon", 0.01);
        Value<Double> t;
        Value<Double> v;
        Value<Double[]> s;
        Value<Number> meanT = new Value<>("menaT", 2.3);
        Value<Number> sdT = new Value<>("sdT", 0.1);
        Value<Number> meanV = new Value<>("meanV", 0.1);
        Value<Number> sdV = new Value<>("sdV", 0.05);
        Value<Number> meanS = new Value<>("meanS", 0.04);
        Value<Number> sdS = new Value<>("sdS", 0.001);
        LogNormal logNormalT = new LogNormal(meanT, sdT, null);
        LogNormal logNormalV = new LogNormal(meanV, sdV,null);
        LogNormal logNormalS = new LogNormal(meanS, sdS, null);


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

        LogNormal logNormal = new LogNormal(M, S, null);
        theta = logNormal.sample();
        Coalescent coalescent = new Coalescent(theta, n, null);
        psi = coalescent.sample();
        PhyloCTMC phyloCTMC = new PhyloCTMC(psi, null, null, Q,
                null, null, l, null, null);
        data = phyloCTMC.sample(); // alignment genotypes


        PloidyModel ploidyModel = new PloidyModel(l= l, n= n, delta = delta);
        Value<Integer[][]> alpha = ploidyModel.sample();

        t = logNormalT.sample();
        v = logNormalV.sample();
        Double[] s1 = new Double[n.value()];
        for (int i = 0; i < n.value(); i++) {
            s1[i] = logNormalS.sample().value();
        }
        s = new Value<>("s", s1);

        CoverageModel coverageModel = new CoverageModel(alpha = alpha, t = t, v= v, s = s);
        Value<Integer[][]> coverage = coverageModel.sample();
        ReadCountModel readCountModel= new ReadCountModel(data, coverage, alpha, epsilon);
        Value<ReadCountData> readCount = readCountModel.sample();



    }

}
