package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulation;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulationFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.GompertzGrowth;


public class GompertzToBEAST implements ValueToBEAST<GompertzPopulation,GompertzGrowth> {

    public GompertzGrowth valueToBEAST(Value<GompertzPopulation> lphyPopFuncVal, BEASTContext context) {

        GompertzGrowth beastPopFunc;

        GompertzPopulationFunction gen = (GompertzPopulationFunction) lphyPopFuncVal.getGenerator();

        RealParameter f0Param = context.getAsRealParameter(gen.getF0());
        RealParameter bParam = context.getAsRealParameter(gen.getB());
        RealParameter NInfinityParam = context.getAsRealParameter(gen.getNInfinity());

        beastPopFunc = new GompertzGrowth();

        beastPopFunc.setInputValue("f0", f0Param);
        beastPopFunc.setInputValue("b", bParam);
        beastPopFunc.setInputValue("NInfinity", NInfinityParam);
        beastPopFunc.initAndValidate();

        //addUpDownOperator(f0Param, bParam, context);

        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        return beastPopFunc;
    }


//    private void addUpDownOperator(RealParameter f0, RealParameter b, BEASTContext context) {
//        UpDownOperator upDownOperator = new UpDownOperator();
//        upDownOperator.setInputValue("scaleFactor", 0.75);
//        upDownOperator.setInputValue("weight", 3.0);
//
//        Tree tree = (Tree) context.getBEASTObject("tree");
//
//        upDownOperator.setInputValue("up", Arrays.asList(f0, b));
//        upDownOperator.setInputValue("down", Arrays.asList(tree));
//        upDownOperator.initAndValidate();
//
//        context.addExtraOperator(upDownOperator);
//    }


    public Class getValueClass() {
        return GompertzPopulation.class;
    }

    public Class<GompertzGrowth> getBEASTClass() {
        return GompertzGrowth.class;
    }

}
