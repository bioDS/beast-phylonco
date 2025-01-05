package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.LogisticPopulation;
import lphy.base.evolution.coalescent.populationmodel.LogisticPopulationFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.LogisticGrowth;

/**
 * Bridges LPhy's LogisticPopulation to BEAST's LogisticGrowth model.
 */
public class LogisticToBEAST implements ValueToBEAST<LogisticPopulation, LogisticGrowth> {

    private static final int I_NA_LOWER_BOUND = 0;
    private static final int I_NA_UPPER_BOUND = 1;

    @Override
    public LogisticGrowth valueToBEAST(Value<LogisticPopulation> lphyPopFuncVal, BEASTContext context) {

        // 1. Get the LogisticPopulationFunction generator from the LPhy Value
        LogisticPopulationFunction gen = (LogisticPopulationFunction) lphyPopFuncVal.getGenerator();

        // 2. Convert LPhy parameters to BEAST RealParameters
        RealParameter bParam = context.getAsRealParameter(gen.getB());
        RealParameter carryingCapacityParam = context.getAsRealParameter(gen.getNCarryingCapacity());
        RealParameter t50Param = context.getAsRealParameter(gen.getT50());

        // 3. Create a new LogisticGrowth model in BEAST
        LogisticGrowth beastPopFunc = new LogisticGrowth();

        // 4. Set required inputs: b, nCarryingCapacity, t50
        beastPopFunc.setInputValue("b", bParam);
        beastPopFunc.setInputValue("nCarryingCapacity", carryingCapacityParam);
        beastPopFunc.setInputValue("t50", t50Param);

        // 5. If NA > 0, convert and set it
        Value<Double> naValue = gen.getNA();
        if (naValue != null && naValue.value() != null && naValue.value() > 0.0) {
            RealParameter naParam = context.getAsRealParameter(naValue);
            beastPopFunc.setInputValue("NA", naParam);
        }


        Value<Integer> iNaValue = gen.getI_na();
        if (iNaValue != null && iNaValue.value() != null) {
            IntegerParameter iNaParam = context.getAsIntegerParameter(iNaValue);
            setIntegerParamBounds(iNaParam, I_NA_LOWER_BOUND, I_NA_UPPER_BOUND);
            beastPopFunc.setInputValue("I_na", iNaParam);
        }


        // 7. Validate and set ID
        beastPopFunc.initAndValidate();
        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        return beastPopFunc;
    }

    /**
     * A helper method to clamp an IntegerParameter between lowerBound and upperBound.
     */
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
    public Class<LogisticPopulation> getValueClass() {
        return LogisticPopulation.class;
    }

    @Override
    public Class<LogisticGrowth> getBEASTClass() {
        return LogisticGrowth.class;
    }
}
