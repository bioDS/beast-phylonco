//package phylonco.lphybeast.tobeast.values;
//
//import beast.base.inference.parameter.RealParameter;
//import lphy.base.evolution.coalescent.PopulationFunction;
//import lphy.base.evolution.coalescent.populationmodel.GompertzPopulation;
//import lphy.base.evolution.coalescent.populationmodel.GompertzPopulationFunction;
//import lphy.core.model.Value;
//import lphybeast.BEASTContext;
//import lphybeast.ValueToBEAST;
//import lphybeast.tobeast.values.ValueToParameter;
//import phylonco.beast.evolution.populationmodel.GompertzGrowth;
//
//public class PopulationFunctionToBEAST implements ValueToBEAST<PopulationFunction, beast.base.evolution.tree.coalescent.PopulationFunction.Abstract> {
//
//    public beast.base.evolution.tree.coalescent.PopulationFunction.Abstract valueToBEAST(Value<PopulationFunction> lphyPopFuncVal, BEASTContext context) {
//        beast.base.evolution.tree.coalescent.PopulationFunction.Abstract beastPopFunc;
//
//        if (lphyPopFuncVal.getType().isAssignableFrom(GompertzPopulation.class)) {
//
//            GompertzPopulationFunction gen = (GompertzPopulationFunction) lphyPopFuncVal.getGenerator();
//
//            RealParameter f0Param = context.getAsRealParameter(gen.getF0());
//            RealParameter bParam = context.getAsRealParameter(gen.getB());
//            RealParameter NInfinityParam = context.getAsRealParameter(gen.getNInfinity());
//
//            beastPopFunc = new GompertzGrowth();
//
//            beastPopFunc.setInputValue("f0", f0Param);
//            beastPopFunc.setInputValue("b", bParam);
//            beastPopFunc.setInputValue("NInfinity", NInfinityParam);
//            beastPopFunc.initAndValidate();
//
//            ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);
//
//            return beastPopFunc;
//        }
//
//        throw new UnsupportedOperationException("TODO");
//    }
//
//    public Class getValueClass() {
//        return PopulationFunction.class;
//    }
//
//    public Class<beast.base.evolution.tree.coalescent.PopulationFunction.Abstract> getBEASTClass() {
//        return beast.base.evolution.tree.coalescent.PopulationFunction.Abstract.class;
//    }
//
//}
