package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import lphy.base.evolution.coalescent.PopulationFunction;
import lphy.base.evolution.coalescent.PopulationFunctionCoalescent;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulation;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import phylonco.beast.evolution.populationmodel.GompertzGrowth;

public class PopFuncCoalescentToBEAST implements
        GeneratorToBEAST<PopulationFunctionCoalescent, beast.base.evolution.tree.coalescent.Coalescent> {
    @Override
    public beast.base.evolution.tree.coalescent.Coalescent generatorToBEAST(PopulationFunctionCoalescent coalescent, BEASTInterface value, BEASTContext context) {

        beast.base.evolution.tree.coalescent.Coalescent beastCoalescent = new beast.base.evolution.tree.coalescent.Coalescent();

        TreeIntervals treeIntervals = new TreeIntervals();
        treeIntervals.setInputValue("tree", value);
        treeIntervals.initAndValidate();

        beastCoalescent.setInputValue("treeIntervals", treeIntervals);

        beast.base.evolution.tree.coalescent.PopulationFunction.Abstract populationFunction;
        //TODO why private ?
        Value<PopulationFunction> lphyPF = coalescent.getParams().get("popFunc");

        if (lphyPF.getType().isAssignableFrom(GompertzPopulation.class)) {

            populationFunction = (GompertzGrowth)  context.getBEASTObject(lphyPF);

        } else {
            // TODO other pop function types
            populationFunction = new ConstantPopulation();
            populationFunction.setInputValue("popSize", context.getBEASTObject("TODO"));
            populationFunction.initAndValidate();
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
