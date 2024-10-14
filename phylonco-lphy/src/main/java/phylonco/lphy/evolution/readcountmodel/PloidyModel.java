package phylonco.lphy.evolution.readcountmodel;

import lphy.base.distribution.Binomial;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.Map;

public class PloidyModel implements GenerativeDistribution<Integer2DMatrix>  {
    private Value<Integer> l;
    private Value<Integer> n;
    private Value<Double> delta;



    public static final String lParamName = "l";
    public static final String nParamName = "n";
    public static final String deltaParamName = "delta";

    public PloidyModel(
            @ParameterInfo(name = lParamName, description = "the length of sequencing.") Value<Integer> l,
            @ParameterInfo(name = nParamName, description = "the number of cells.") Value<Integer> n,
            @ParameterInfo(name = deltaParamName, description = "allelic dropout error probability.") Value<Double> delta

            ){
        super();
        this.l = l;
        this.n = n;
        this.delta = delta;

    }

    @GeneratorInfo(
            name = "Ploidy",
            narrativeName = "ploidy model",
            description = "Observed ploidy after allelic dropout.")


    @Override
    public RandomVariable<Integer2DMatrix> sample() {
        RandomVariable<Integer2DMatrix> alpha;
        Value<Integer> numberA = new Value<>("numberA", 1);
        Binomial binomialAlpha = new Binomial(delta, numberA);

        Integer[][] alp = new Integer[this.n.value()][this.l.value()];
        for (int i = 0; i < this.n.value(); i++) {
            for (int j = 0; j < this.l.value(); j++) {
                Value<Integer> alpha1 = binomialAlpha.sample();
                if (alpha1.value() == 1) {
                    alp[i][j] = 1;
                } else alp[i][j] = 2;
            }
        }

        alpha = new RandomVariable<>("alpha", new Integer2DMatrix(alp), this);
        return alpha;
    }


    @Override
    public Map<String, Value> getParams() {
        return Map.of(
                lParamName,l,
                nParamName,n,
                deltaParamName,delta
        );
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(lParamName)) l = value;
        else if (paramName.equals(nParamName)) n = value;
        else if (paramName.equals(deltaParamName)) delta = value;
        else throw new RuntimeException("Unrecognised parameter name: " + paramName);
    }

    public Value<Integer> getL() {
        return l;
    }

    public Value<Integer> getN() {
        return n;
    }

    public Value<Double> getDelta() {
        return delta;
    }
}
