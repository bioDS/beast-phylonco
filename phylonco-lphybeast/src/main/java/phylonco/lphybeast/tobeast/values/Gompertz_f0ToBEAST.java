package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulationFunction_f0;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulation_f0;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.GompertzGrowth_f0;

/**
 * Bridges LPhy's GompertzPopulation_f0 to BEAST's GompertzGrowth_f0,
 * handling optional NA via I_na just like ExponentialToBEAST.
 */
public class Gompertz_f0ToBEAST implements ValueToBEAST<GompertzPopulation_f0, GompertzGrowth_f0> {

    private static final int I_NA_LOWER_BOUND = 0;
    private static final int I_NA_UPPER_BOUND = 1;

    @Override
    public GompertzGrowth_f0 valueToBEAST(Value<GompertzPopulation_f0> lphyPopFuncVal, BEASTContext context) {

        GompertzPopulationFunction_f0 gen = (GompertzPopulationFunction_f0) lphyPopFuncVal.getGenerator();

        RealParameter f0Param = context.getAsRealParameter(gen.getF0());
        RealParameter bParam  = context.getAsRealParameter(gen.getB());
        RealParameter N0Param = context.getAsRealParameter(gen.getN0());

        GompertzGrowth_f0 beastPopFunc = new GompertzGrowth_f0();
        beastPopFunc.setInputValue("f0", f0Param);
        beastPopFunc.setInputValue("b", bParam);
        beastPopFunc.setInputValue("N0", N0Param);

        Value<Double> naValue = gen.getNA();
        if (naValue != null && naValue.value() != null && naValue.value() > 0.0) {
            RealParameter NAParam = context.getAsRealParameter(naValue);
            beastPopFunc.setInputValue("NA", NAParam);
        }

        Value<Integer> iNaValue = gen.getI_na();
        if (iNaValue != null && iNaValue.value() != null) {
            IntegerParameter iNaParam = context.getAsIntegerParameter(iNaValue);
            setIntegerParamBounds(iNaParam, I_NA_LOWER_BOUND, I_NA_UPPER_BOUND);
            beastPopFunc.setInputValue("I_na", iNaParam);
        }

        beastPopFunc.initAndValidate();
        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        return beastPopFunc;
    }

    private void setIntegerParamBounds(IntegerParameter param, int lowerBound, int upperBound) {
        param.setInputValue("lower", lowerBound);
        param.setInputValue("upper", upperBound);

        int currentValue = param.getValue();
        if (currentValue < lowerBound) {
            param.setValue(lowerBound);
        } else if (currentValue > upperBound) {
            param.setValue(upperBound);
        }
    }

    @Override
    public Class<GompertzPopulation_f0> getValueClass() {
        return GompertzPopulation_f0.class;
    }

    @Override
    public Class<GompertzGrowth_f0> getBEASTClass() {
        return GompertzGrowth_f0.class;
    }
}
