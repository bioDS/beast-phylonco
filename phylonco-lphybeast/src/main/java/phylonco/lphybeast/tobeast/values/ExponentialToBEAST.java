package phylonco.lphybeast.tobeast.values;

import beast.base.evolution.tree.coalescent.ExponentialGrowth;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.ExponentialPopulation;
import lphy.base.evolution.coalescent.populationmodel.ExponentialPopulationFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;

public class ExponentialToBEAST implements ValueToBEAST<ExponentialPopulation, ExponentialGrowth> {

    public ExponentialGrowth valueToBEAST(Value<ExponentialPopulation> lphyPopFuncVal, BEASTContext context) {

        ExponentialGrowth beastPopFunc;

        ExponentialPopulationFunction gen = (ExponentialPopulationFunction) lphyPopFuncVal.getGenerator();

        RealParameter GrowthRateParam = context.getAsRealParameter(gen.getGrowthRate());
        RealParameter N0Param = context.getAsRealParameter(gen.getN0());

        beastPopFunc = new ExponentialGrowth();

        beastPopFunc.setInputValue("GrowthRate", GrowthRateParam);
        beastPopFunc.setInputValue("N0", N0Param);
        beastPopFunc.initAndValidate();

        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        return beastPopFunc;
    }

    public Class getValueClass() {
        return ExponentialPopulation.class;
    }

    public Class<ExponentialGrowth> getBEASTClass() {
        return ExponentialGrowth.class;
    }

}