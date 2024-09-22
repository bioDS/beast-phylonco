package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.ExpansionPopulation;
import lphy.base.evolution.coalescent.populationmodel.ExpansionPopulationFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.ExpansionGrowth;


public class ExpansionToBEAST implements ValueToBEAST<ExpansionPopulation, ExpansionGrowth> {

    public ExpansionGrowth valueToBEAST(Value<ExpansionPopulation> lphyPopFuncVal, BEASTContext context) {

        ExpansionGrowth beastPopFunc;

        ExpansionPopulationFunction gen = (ExpansionPopulationFunction) lphyPopFuncVal.getGenerator();

        RealParameter TauParam = context.getAsRealParameter(gen.getTau());
        RealParameter RParam = context.getAsRealParameter(gen.getR());
        RealParameter NCParam = context.getAsRealParameter(gen.getNC());
        RealParameter N0Param = context.getAsRealParameter(gen.getN0());

        beastPopFunc = new ExpansionGrowth();

        beastPopFunc.setInputValue("r", RParam);
        beastPopFunc.setInputValue("N0", N0Param);
        beastPopFunc.setInputValue("NC", NCParam);
        beastPopFunc.setInputValue("tau", TauParam);


        beastPopFunc.initAndValidate();

        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        return beastPopFunc;
    }

    public Class getValueClass() {
        return ExpansionPopulation.class;
    }

    public Class<ExpansionGrowth> getBEASTClass() {
        return ExpansionGrowth.class;
    }

}