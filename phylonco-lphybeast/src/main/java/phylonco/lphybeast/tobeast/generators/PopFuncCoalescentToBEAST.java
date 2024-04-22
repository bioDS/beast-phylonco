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
            UpDownOperator upDownOperator = popFuncWithUpDownOp.getUpDownOperator((Tree) value);
            context.addExtraOperator(upDownOperator);
        }

        beastCoalescent.setInputValue("populationModel", populationFunction);

        beastCoalescent.initAndValidate();

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
