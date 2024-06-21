package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulationFunction_f0;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulation_f0;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.GompertzGrowth_f0;

public class Gompertz_f0ToBEAST implements ValueToBEAST<GompertzPopulation_f0, GompertzGrowth_f0> {

    public GompertzGrowth_f0 valueToBEAST(Value<GompertzPopulation_f0> lphyPopFuncVal, BEASTContext context) {

        GompertzGrowth_f0 beastPopFunc;

        GompertzPopulationFunction_f0 gen = (GompertzPopulationFunction_f0) lphyPopFuncVal.getGenerator();

        RealParameter f0Param = context.getAsRealParameter(gen.getF0());

//        RealParameter N0Param = context.getAsRealParameter(gen.getN0());
        RealParameter bParam = context.getAsRealParameter(gen.getB());
 //       RealParameter NInfinityParam = context.getAsRealParameter(gen.getNInfinity());
        RealParameter N0Param = context.getAsRealParameter(gen.getN0());
//        RealParameter alphaParam = context.getAsRealParameter(gen.getAlpha());
//        RealParameter betaParam = context.getAsRealParameter(gen.getBeta());

//        RealParameter alphaParam = context.getAsRealParameter(gen.getAlpha());
//        RealParameter betaParam = context.getAsRealParameter(gen.getBeta());



        beastPopFunc = new GompertzGrowth_f0();


//        beastPopFunc.setInputValue("N0", N0Param);

        beastPopFunc.setInputValue("f0", f0Param);
        beastPopFunc.setInputValue("b", bParam);
        beastPopFunc.setInputValue("N0", N0Param);
        beastPopFunc.initAndValidate();

        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);


        return beastPopFunc;
    }




    public Class getValueClass() {
        return GompertzPopulation_f0.class;
    }

    public Class<GompertzGrowth_f0> getBEASTClass() {
        return GompertzGrowth_f0.class;
    }

}
