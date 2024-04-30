package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.LogisticPopulation;
import lphy.base.evolution.coalescent.populationmodel.LogisticPopulationFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.LogisticGrowth;

public class LogisticToBEAST implements ValueToBEAST <LogisticPopulation, LogisticGrowth> {
    public LogisticGrowth valueToBEAST(Value<LogisticPopulation> lphyPopFuncVal, BEASTContext context) {

        LogisticGrowth beastPopFunc;

        LogisticPopulationFunction gen = (LogisticPopulationFunction) lphyPopFuncVal.getGenerator();

        RealParameter bParam = context.getAsRealParameter(gen.getB());
        RealParameter NCarryingCapacityParam = context.getAsRealParameter(gen.getNCarryingCapacity());
        RealParameter T50Param = context.getAsRealParameter(gen.getT50());

        beastPopFunc = new LogisticGrowth();

        beastPopFunc.setInputValue("nCarryingCapacity", NCarryingCapacityParam);
        beastPopFunc.setInputValue("b", bParam);
        beastPopFunc.setInputValue("t50", T50Param);
        beastPopFunc.initAndValidate();

        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        return beastPopFunc;

    }

    public Class getValueClass() {
        return LogisticPopulation.class;
    }

    public Class<LogisticGrowth> getBEASTClass() {
        return LogisticGrowth.class;
    }

}
