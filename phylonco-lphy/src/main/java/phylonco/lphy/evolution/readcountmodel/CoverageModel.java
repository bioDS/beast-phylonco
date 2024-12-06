package phylonco.lphy.evolution.readcountmodel;

import lphy.base.distribution.GeneralNegativeBinomial;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.Map;

public class CoverageModel implements GenerativeDistribution<Integer2DMatrix> {
    private Value<Integer2DMatrix> alpha;
    private Value<Double> t;
    private Value<Double> v;
    private Value<Double[]> s;





    public static final String alphaParamName = "alpha";
    public static final String tParamName = "t";
    public static final String vParamName = "v";
    public static final String sParamName = "s";


    public CoverageModel(
            @ParameterInfo(name = alphaParamName, description = "alpha, allelic dropout events for each cell at each site.") Value<Integer2DMatrix> alpha,
            @ParameterInfo(name = tParamName, description = "t, mean of allelic coverage.") Value<Double> t,
            @ParameterInfo(name = vParamName, description = "v, variance of allelic coverage.") Value<Double> v,
            @ParameterInfo(name = sParamName, description = "s, size factor of cell.") Value<Double[]> s

            ) {
        super();
        this.alpha = alpha;
        this.t = t;
        this.v = v;
        this.s = s;


    }


    @GeneratorInfo(
            name = "CoverageModel",
            narrativeName = "coverage model",
            description = "A model to simulate the coverage at each site in each cell")
    @Override
    public RandomVariable<Integer2DMatrix> sample(){
        RandomVariable<Integer2DMatrix> coverage;
        int n = alpha.value().nTaxa();
        int l = alpha.value().nchar();
        Value<Double> r;
        Value<Double> p;
        Integer[][] cov = new Integer[n][l];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < l; j++) {
                Double mean = (this.alpha.value().getState(i,j) * this.t.value() * this.s.value()[i]);
                Double variance = mean + (Math.pow((double) alpha.value().getState(i,j), 2.0) * v.value() * Math.pow(this.s.value()[i], 2.0));
                double pValue = mean / variance;
                double rValue = Math.pow(mean, 2) / (variance - mean);
                r = new Value<>("r", rValue);
                p = new Value<>("p", pValue);
                if (rValue == 0){
                    cov[i][j] = 0;
                }else {
                    GeneralNegativeBinomial negativeBinomial = new GeneralNegativeBinomial(r, p);
                    cov[i][j] = negativeBinomial.sample().value();
                }
            }
        }
        coverage = new RandomVariable<>("coverage", new Integer2DMatrix(cov), this);
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

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(alphaParamName)) alpha = value;
        else if (paramName.equals(tParamName)) t = value;
        else if (paramName.equals(vParamName)) v = value;
        else if (paramName.equals(sParamName)) s = value;
        else throw new RuntimeException("Unrecognised parameter name: " + paramName);
    }

}
