package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.ConstantPopulation;
import lphy.base.evolution.coalescent.populationmodel.ConstantPopulationFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.ConstantGrowth;


public class ConstantToBEAST implements ValueToBEAST<ConstantPopulation, ConstantGrowth> {

    public ConstantGrowth valueToBEAST(Value<ConstantPopulation> lphyPopFuncVal, BEASTContext context) {

        ConstantGrowth beastPopFunc;

        ConstantPopulationFunction gen = (ConstantPopulationFunction) lphyPopFuncVal.getGenerator();


        RealParameter N0Param = context.getAsRealParameter(gen.getN0());

        beastPopFunc = new ConstantGrowth();


        beastPopFunc.setInputValue("N0", N0Param);
        beastPopFunc.initAndValidate();

        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        return beastPopFunc;
    }

    public Class getValueClass() {
        return ConstantPopulation.class;
    }

    public Class<ConstantGrowth> getBEASTClass() {
        return ConstantGrowth.class;
    }

}