package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.ExponentialPopulation;
import lphy.base.evolution.coalescent.populationmodel.ExponentialPopulationFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.ExponentialGrowth;

/**
 * Bridges LPhy's ExponentialPopulation to BEAST's ExponentialGrowth model.
 */
public class ExponentialToBEAST implements ValueToBEAST<ExponentialPopulation, ExponentialGrowth> {

    // Optionally define constants for bounds
    private static final int I_NA_LOWER_BOUND = 0;
    private static final int I_NA_UPPER_BOUND = 1;

    @Override
    public ExponentialGrowth valueToBEAST(Value<ExponentialPopulation> lphyPopFuncVal, BEASTContext context) {

        ExponentialPopulationFunction gen = (ExponentialPopulationFunction) lphyPopFuncVal.getGenerator();

        RealParameter growthRateParam = context.getAsRealParameter(gen.getGrowthRate());
        RealParameter n0Param = context.getAsRealParameter(gen.getN0());

        ExponentialGrowth beastPopFunc = new ExponentialGrowth();

        beastPopFunc.setInputValue("GrowthRate", growthRateParam);
        beastPopFunc.setInputValue("N0", n0Param);

        Value<Double> naValue = gen.getNA();
        if (naValue != null && naValue.value() != null && naValue.value() > 0.0) {
            RealParameter NAParam = context.getAsRealParameter(naValue);
            beastPopFunc.setInputValue("NA", NAParam);
        }

        Value<Integer> iNaValue = gen.getI_na();
        IntegerParameter iNaParam = context.getAsIntegerParameter(iNaValue);
        setIntegerParamBounds(iNaParam, I_NA_LOWER_BOUND, I_NA_UPPER_BOUND);

        beastPopFunc.setInputValue("I_na", iNaParam);

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
    public Class<ExponentialPopulation> getValueClass() {
        return ExponentialPopulation.class;
    }

    @Override
    public Class<ExponentialGrowth> getBEASTClass() {
        return ExponentialGrowth.class;
    }
}
