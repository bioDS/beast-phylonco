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

        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);


        return beastPopFunc;
    }

//    private void addScaleOperators(RealParameter f0, RealParameter b, RealParameter NInfinity, BEASTContext context) {
//        ScaleOperator f0Scaler = new ScaleOperator();
//        f0Scaler.setInputValue("parameter", f0);
//        f0Scaler.setInputValue("scaleFactor", 0.75);
//        f0Scaler.setInputValue("weight", 0.1);
//        f0Scaler.initAndValidate();
//
//        ScaleOperator bScaler = new ScaleOperator();
//        bScaler.setInputValue("parameter", b);
//        bScaler.setInputValue("scaleFactor", 0.75);
//        bScaler.setInputValue("weight", 0.1);
//        bScaler.initAndValidate();
//
//        ScaleOperator NInfinityScaler = new ScaleOperator();
//        NInfinityScaler.setInputValue("parameter", NInfinity);
//        NInfinityScaler.setInputValue("scaleFactor", 0.75);
//        NInfinityScaler.setInputValue("weight", 0.1);
//        NInfinityScaler.initAndValidate();
//
//        // Assuming context is a facility to add these operators to the model
//        context.addExtraOperator(f0Scaler);
//        context.addExtraOperator(bScaler);
//        context.addExtraOperator(NInfinityScaler);
//    }



    public Class getValueClass() {
        return GompertzPopulation.class;
    }

    public Class<GompertzGrowth> getBEASTClass() {
        return GompertzGrowth.class;
    }

}
