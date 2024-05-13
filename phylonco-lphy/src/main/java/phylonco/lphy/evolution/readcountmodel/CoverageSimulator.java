package phylonco.lphy.evolution.readcountmodel;

import lphy.base.distribution.NegativeBinomial;
import lphy.base.evolution.alignment.Alignment;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.ParameterInfo;

import java.util.Map;

public class CoverageSimulator implements GenerativeDistribution<Integer[][]> {
    private NegativeBinomial negativeBinomial;
    private Value<Integer[][]> alpha;
    private Value<Integer> t;
    private Value<Integer> v;
    private Value<Integer[]> s;
    private RandomVariable<Integer[][]> coverage;




    public static final String alphaParamName = "alpha";
    public static final String tParamName = "t";
    public static final String vParamName = "v";
    public static final String sParamName = "s";


    public CoverageSimulator(
            @ParameterInfo(name = alphaParamName, description = "alpha, allelic dropout events for each cell at each site.") Value<Integer[][]> alpha,
            @ParameterInfo(name = tParamName, description = "t") Value<Integer> t,
            @ParameterInfo(name = vParamName, description = "v") Value<Integer> v,
            @ParameterInfo(name = sParamName, description = "s") Value<Integer[]> s

            ) {
        super();
        this.alpha = alpha;
        this.t = t;
        this.v = v;
        this.s = s;


    }

    @Override
    public RandomVariable<Integer[][]> sample(){
        int n = alpha.value().length;
        int l = alpha.value()[0].length;
        Value<Integer> r;
        Value<Double> p;
        Integer[][] cov = new Integer[n][l];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < l; j++) {
                Double mean = (double) (this.alpha.value()[i][j] * this.t.value() * this.s.value()[i]);
                Double variance = mean + (Math.pow((double) alpha.value()[i][j], 2.0) * v.value() * Math.pow(this.s.value()[i], 2.0));
                double pValue = mean / variance;
                float rFloat = (float) (Math.pow(mean, 2) / (variance - mean));
                int rValue = Math.round(rFloat);
                r = new Value<>("r", rValue);
                p = new Value<>("p", pValue);
                NegativeBinomial negativeBinomial = new NegativeBinomial(r, p);
                cov[i][j] = negativeBinomial.sample().value();
            }
        }
        coverage = new RandomVariable<>("coverage", cov, this);
        return coverage;
    }

    @Override
    public Map<String, Value> getParams() {
        return Map.of(
                alphaParamName, alpha,
                tParamName, t,
                vParamName, v,
                sParamName, s
        );
    }

}
