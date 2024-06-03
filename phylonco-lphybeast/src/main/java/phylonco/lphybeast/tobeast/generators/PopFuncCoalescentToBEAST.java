package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.inference.operator.UpDownOperator;
import lphy.base.evolution.coalescent.PopulationFunction;
import lphy.base.evolution.coalescent.PopulationFunctionCoalescent;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import phylonco.beast.evolution.populationmodel.PopFuncWithUpDownOp;

public class PopFuncCoalescentToBEAST implements
        GeneratorToBEAST<PopulationFunctionCoalescent, beast.base.evolution.tree.coalescent.Coalescent> {
    @Override
    public beast.base.evolution.tree.coalescent.Coalescent generatorToBEAST(PopulationFunctionCoalescent coalescent, BEASTInterface value, BEASTContext context) {

        beast.base.evolution.tree.coalescent.Coalescent beastCoalescent = new beast.base.evolution.tree.coalescent.Coalescent();

        Value<PopulationFunction> lphyPF = coalescent.getParams().get("popFunc");
        if (lphyPF == null) {
            throw new IllegalArgumentException("popFunc parameter is null.");
        }

        beast.base.evolution.tree.coalescent.PopulationFunction populationFunction =
                (beast.base.evolution.tree.coalescent.PopulationFunction) context.getBEASTObject(lphyPF);

        TreeIntervals treeIntervals = new TreeIntervals();
        treeIntervals.setInputValue("tree", value);
        treeIntervals.initAndValidate();



        beastCoalescent.setInputValue("treeIntervals", treeIntervals);

        //        beast.base.evolution.tree.coalescent.PopulationFunction.Abstract populationFunction;
        //
        //        if (lphyPF.getType().isAssignableFrom(GompertzPopulation.class)) {
        //
        //            populationFunction = (GompertzGrowth)  context.getBEASTObject(lphyPF);
        //
        //        } else {
        //            // TODO other pop function types
        //            populationFunction = new ConstantPopulation();
        //            populationFunction.setInputValue("popSize", context.getBEASTObject("TODO"));
        //            populationFunction.initAndValidate();
        //        }

        if (populationFunction instanceof PopFuncWithUpDownOp popFuncWithUpDownOp) {
            UpDownOperator upDownOperator1 = popFuncWithUpDownOp.getUpDownOperator1((Tree) value);
            UpDownOperator upDownOperator2 = popFuncWithUpDownOp.getUpDownOperator2((Tree) value);
            context.addExtraOperator(upDownOperator1);
            context.addExtraOperator(upDownOperator2);
        }

        beastCoalescent.setInputValue("populationModel", populationFunction);

        beastCoalescent.initAndValidate();


//        AdaptableVarianceMultivariateNormalOperator avmnOperator = new AdaptableVarianceMultivariateNormalOperator();
        ////        avmnOperator.setID("AVMNOperator");
        ////        avmnOperator.setInputValue("beta", 0.05);
        ////        avmnOperator.setInputValue("burnin", 400);
        ////        avmnOperator.setInputValue("initial", 800);
        ////        avmnOperator.setInputValue("weight", 2.0);
        ////
        ////
        ////        // Set up transformations
        ////        List<Transform> transformations = new ArrayList<>();
        ////
        ////        // f0 transformation
        ////        RealParameter f0Param = (RealParameter) context.getBEASTObject((Value<RealParameter>) coalescent.getParams().get("f0"));
        ////        if (f0Param == null) {
        ////            throw new IllegalArgumentException("f0 parameter is null.");
        ////        }
        ////        Transform.Interval f0Transform = new Transform.Interval();
        ////        f0Transform.setInputValue("lower", 0.0);
        ////        f0Transform.setInputValue("upper", 1.0);
        ////        f0Transform.setInputValue("f", f0Param);
        ////        transformations.add(f0Transform);
        ////
        ////        RealParameter bParam = (RealParameter) context.getBEASTObject((Value<RealParameter>) coalescent.getParams().get("b"));
        ////        if (bParam == null) {
        ////            throw new IllegalArgumentException("b parameter is null.");
        ////        }
        ////        Transform.LogTransform bTransform = new Transform.LogTransform();
        ////        bTransform.setInputValue("f", bParam);
        ////        transformations.add(bTransform);
        ////
        ////        avmnOperator.setInputValue("transformations", transformations);
        ////        avmnOperator.initAndValidate();
        ////
        ////        context.addExtraOperator(avmnOperator);

        return beastCoalescent;
    }







    @Override
    public Class<PopulationFunctionCoalescent> getGeneratorClass() {
        return PopulationFunctionCoalescent.class;
    }

    @Override
    public Class<beast.base.evolution.tree.coalescent.Coalescent> getBEASTClass() {
        return beast.base.evolution.tree.coalescent.Coalescent.class;
    }
}