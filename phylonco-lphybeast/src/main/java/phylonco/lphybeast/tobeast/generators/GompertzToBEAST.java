//package phylonco.lphybeast.tobeast.generators;
//
//import beast.base.core.BEASTInterface;
//import beast.base.core.Function;
//import beast.base.evolution.tree.coalescent.PopulationFunction;
//import beast.base.inference.parameter.RealParameter;
//import lphy.base.evolution.coalescent.populationmodel.GompertzPopulation;
//import lphy.base.evolution.coalescent.populationmodel.GompertzPopulationFunction;
//import lphy.core.model.Generator;
//import lphy.core.model.ValueUtils;
//import lphybeast.BEASTContext;
//import lphybeast.GeneratorToBEAST;
//import phylonco.beast.evolution.populationmodel.GompertzGrowth;
//import phylonco.lphy.evolution.substitutionmodel.GT16;
//
//import java.util.List;
//
//public class GompertzToBEAST implements GeneratorToBEAST<GompertzPopulationFunction, GompertzGrowth> {
////     <populationModel id="gompertzPopulationModel" spec="phylonco.beast.evolution.populationmodel.GompertzGrowth">
////                        <parameter name="f0" idref="f0"/>
////                        <parameter name="NInfinity" idref="NInfinity"/>
////                        <parameter name="b" idref="b"/>
////                    </populationModel>
//
//    @Override
//    public GompertzGrowth generatorToBEAST(GompertzPopulationFunction gompertzPopulationFunction, BEASTInterface beastInterface, BEASTContext beastContext) {
//
//        GompertzGrowth beastGompertzGrowth = new GompertzGrowth();
//
//        RealParameter f0Param = beastContext.getAsRealParameter(gompertzPopulationFunction.getF0());
//        RealParameter bParam = beastContext.getAsRealParameter(gompertzPopulationFunction.getB());
//        RealParameter NInfinityParam = beastContext.getAsRealParameter(gompertzPopulationFunction.getNInfinity());
//
//        beastGompertzGrowth.setInputValue("f0", f0Param);
//        beastGompertzGrowth.setInputValue("b", bParam);
//        beastGompertzGrowth.setInputValue("NInfinity", NInfinityParam);
//        beastGompertzGrowth.initAndValidate();
//
//        return beastGompertzGrowth;
//    }
//
//    @Override
//    public Class<GompertzPopulationFunction> getGeneratorClass() { return GompertzPopulationFunction.class; }
//
//    @Override
//    public Class<GompertzGrowth> getBEASTClass() {
//        return GompertzGrowth.class;
//    }
//}
