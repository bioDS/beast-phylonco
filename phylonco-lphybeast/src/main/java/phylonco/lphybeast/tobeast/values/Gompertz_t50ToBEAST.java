package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulationFunction_t50;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulation_t50;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.GompertzGrowth_t50;

public class Gompertz_t50ToBEAST implements ValueToBEAST<GompertzPopulation_t50, GompertzGrowth_t50> {

    public GompertzGrowth_t50 valueToBEAST(Value<GompertzPopulation_t50> lphyPopFuncVal, BEASTContext context) {

        GompertzGrowth_t50 beastPopFunc;

        GompertzPopulationFunction_t50 gen = (GompertzPopulationFunction_t50) lphyPopFuncVal.getGenerator();

        RealParameter NInfinityParam = context.getAsRealParameter(gen.getNInfinity());

        //        RealParameter N0Param = context.getAsRealParameter(gen.getN0());
        RealParameter bParam = context.getAsRealParameter(gen.getB());
        //       RealParameter NInfinityParam = context.getAsRealParameter(gen.getNInfinity());
        RealParameter t50Param = context.getAsRealParameter(gen.getT50());
        //        RealParameter alphaParam = context.getAsRealParameter(gen.getAlpha());
        //        RealParameter betaParam = context.getAsRealParameter(gen.getBeta());

        //        RealParameter alphaParam = context.getAsRealParameter(gen.getAlpha());
        //        RealParameter betaParam = context.getAsRealParameter(gen.getBeta());



        beastPopFunc = new GompertzGrowth_t50();


        //        beastPopFunc.setInputValue("N0", N0Param);

        beastPopFunc.setInputValue("t50", t50Param);
        beastPopFunc.setInputValue("b", bParam);
        beastPopFunc.setInputValue("NInfinity", NInfinityParam);
        beastPopFunc.initAndValidate();

        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);


        return beastPopFunc;
    }




    public Class getValueClass() {
        return GompertzPopulation_t50.class;
    }

    public Class<GompertzGrowth_t50> getBEASTClass() {
        return GompertzGrowth_t50.class;
    }

}

