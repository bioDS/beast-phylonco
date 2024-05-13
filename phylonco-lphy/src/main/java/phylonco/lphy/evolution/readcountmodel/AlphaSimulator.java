package phylonco.lphy.evolution.readcountmodel;

import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.ParameterInfo;

import java.util.Map;

public class AlphaSimulator implements GenerativeDistribution<Integer[][]>  {
    private Value<Integer> l;
    private Value<Integer> n;
    private Value<Double> delta;

    private RandomVariable<Integer[][]> alpha;



    public static final String lParamName = "l";
    public static final String nParamName = "n";
    public static final String deltaParamName = "delta";




    private AlphaSimulator(
            @ParameterInfo(name = lParamName, description = "the length of sequencing") Value<Integer> l,
            @ParameterInfo(name = nParamName, description = "the number of cells.") Value<Integer> n,
            @ParameterInfo(name = deltaParamName, description = "allelic dropout error probability.") Value<Double> delta

            ){
        super();
        this.l = l;
        this.n = n;
        this.delta = delta;



    }

    @Override
    public RandomVariable<Integer[][]> sample() {
        Multinomial multinomialAlpha;
        Double[] proAlpha = {this.delta.value(), 1-this.delta.value()};
        Value<Double[]> proA = new Value<>("proA", proAlpha);
        Value<Integer> numberA = new Value<>("numberA", 1);
        multinomialAlpha = new Multinomial(numberA, proA);

        Integer[][] alp = new Integer[n.value()][l.value()];
        for (int i = 0; i < n.value(); i++) {
            for (int j = 0; j < l.value(); j++) {
                Value<Integer[]> alpha1 = multinomialAlpha.sample();
                if (alpha1.value()[0] == 1) {
                    alp[i][j] = 1;
                } else alp[i][j] = 2;
            }
        }

        alpha = new RandomVariable<>("alpha", alp, this);
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



}
