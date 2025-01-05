package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulationFunction_t50;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulation_t50;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.GompertzGrowth_t50;

/**
 * Bridges an LPhy GompertzPopulation_t50 (via GompertzPopulationFunction_t50)
 * to a BEAST GompertzGrowth_t50, including an optional ancestral population size NA
 * and indicator I_na (0 or 1).
 */
public class Gompertz_t50ToBEAST implements ValueToBEAST<GompertzPopulation_t50, GompertzGrowth_t50> {

    private static final int I_NA_LOWER_BOUND = 0;
    private static final int I_NA_UPPER_BOUND = 1;

    @Override
    public GompertzGrowth_t50 valueToBEAST(Value<GompertzPopulation_t50> lphyPopFuncVal, BEASTContext context) {

        GompertzPopulationFunction_t50 gen = (GompertzPopulationFunction_t50) lphyPopFuncVal.getGenerator();

        RealParameter t50Param       = context.getAsRealParameter(gen.getT50());
        RealParameter bParam         = context.getAsRealParameter(gen.getB());
        RealParameter NInfinityParam = context.getAsRealParameter(gen.getNInfinity());

        RealParameter NAParam = null;
        if (gen.getNA() != null && gen.getNA().value() != null) {
            NAParam = context.getAsRealParameter(gen.getNA());
        }

        GompertzGrowth_t50 beastPopFunc = new GompertzGrowth_t50();
        beastPopFunc.setInputValue("t50",       t50Param);
        beastPopFunc.setInputValue("b",         bParam);
        beastPopFunc.setInputValue("NInfinity", NInfinityParam);

        if (NAParam != null) {
            beastPopFunc.setInputValue("NA", NAParam);
        }

        Value<Integer> iNaValue = gen.getI_na();
        if (iNaValue != null && iNaValue.value() != null) {
            IntegerParameter iNaParam = context.getAsIntegerParameter(iNaValue);
            clampIntegerParam(iNaParam, I_NA_LOWER_BOUND, I_NA_UPPER_BOUND);
            beastPopFunc.setInputValue("I_na", iNaParam);
        }

        beastPopFunc.initAndValidate();
        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        return beastPopFunc;
    }

    /**
     * Clamps an IntegerParameter value to [lower, upper].
     */
    private void clampIntegerParam(IntegerParameter param, int lower, int upper) {
        param.setInputValue("lower", lower);
        param.setInputValue("upper", upper);

        int currentValue = param.getValue();
        if (currentValue < lower) {
            param.setValue(lower);
        } else if (currentValue > upper) {
            param.setValue(upper);
        }
    }

    @Override
    public Class<GompertzPopulation_t50> getValueClass() {
        return GompertzPopulation_t50.class;
    }

    @Override
    public Class<GompertzGrowth_t50> getBEASTClass() {
        return GompertzGrowth_t50.class;
    }
}
