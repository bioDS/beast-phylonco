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

//        RealParameter N0Param = context.getAsRealParameter(gen.getN0());
        RealParameter bParam = context.getAsRealParameter(gen.getB());
 //       RealParameter NInfinityParam = context.getAsRealParameter(gen.getNInfinity());
        RealParameter N0Param = context.getAsRealParameter(gen.getN0());
//        RealParameter alphaParam = context.getAsRealParameter(gen.getAlpha());
//        RealParameter betaParam = context.getAsRealParameter(gen.getBeta());

//        RealParameter alphaParam = context.getAsRealParameter(gen.getAlpha());
//        RealParameter betaParam = context.getAsRealParameter(gen.getBeta());



        beastPopFunc = new GompertzGrowth();


//        beastPopFunc.setInputValue("N0", N0Param);

        beastPopFunc.setInputValue("f0", f0Param);
        beastPopFunc.setInputValue("b", bParam);
        beastPopFunc.setInputValue("N0", N0Param);
        beastPopFunc.initAndValidate();

        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);


        return beastPopFunc;
    }




    public Class getValueClass() {
        return GompertzPopulation.class;
    }

    public Class<GompertzGrowth> getBEASTClass() {
        return GompertzGrowth.class;
    }

}
