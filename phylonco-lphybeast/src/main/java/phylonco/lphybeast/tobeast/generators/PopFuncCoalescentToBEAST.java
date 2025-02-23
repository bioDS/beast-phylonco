package phylonco.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.operator.kernel.AdaptableVarianceMultivariateNormalOperator;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.inference.operator.UpDownOperator;
import lphy.base.evolution.coalescent.PopulationFunction;
import lphy.base.evolution.coalescent.PopulationFunctionCoalescent;
import lphy.base.evolution.coalescent.populationmodel.SVSPopulation;
import lphy.base.evolution.coalescent.populationmodel.SVSPopulationFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;
import phylonco.beast.evolution.populationmodel.PopFuncWithAVMNOp;
import phylonco.beast.evolution.populationmodel.PopFuncWithUpDownOp;
import phylonco.beast.evolution.populationmodel.PopFuncWithUpOp;

public class PopFuncCoalescentToBEAST implements
        GeneratorToBEAST<PopulationFunctionCoalescent, beast.base.evolution.tree.coalescent.Coalescent> {
    @Override
    public beast.base.evolution.tree.coalescent.Coalescent generatorToBEAST(PopulationFunctionCoalescent coalescent, BEASTInterface value, BEASTContext context) {

        beast.base.evolution.tree.coalescent.Coalescent beastCoalescent = new beast.base.evolution.tree.coalescent.Coalescent();

        Value<PopulationFunction> lphyPF = coalescent.getParams().get("popFunc");
        if (lphyPF == null) {
            throw new IllegalArgumentException("popFunc parameter is null.");
        }

        Object beastObj = context.getBEASTObject(lphyPF);
        if (beastObj == null) {
            throw new IllegalArgumentException("BEAST object conversion returned null for population function.");
        }

        beast.base.evolution.tree.coalescent.PopulationFunction populationFunction;

        if (beastObj instanceof SVSPopulation) {
            populationFunction = convertSVSPopulationFunction((SVSPopulationFunction) beastObj, context);
        } else if (beastObj instanceof beast.base.evolution.tree.coalescent.PopulationFunction) {
            populationFunction = (beast.base.evolution.tree.coalescent.PopulationFunction) beastObj;
        } else {
            throw new IllegalArgumentException("Invalid type for population function: " + beastObj.getClass().getName());
        }

        TreeIntervals treeIntervals = new TreeIntervals();
        treeIntervals.setInputValue("tree", value);
        treeIntervals.initAndValidate();

        beastCoalescent.setInputValue("treeIntervals", treeIntervals);

        if (populationFunction instanceof PopFuncWithUpDownOp) {
            PopFuncWithUpDownOp popFuncWithUpDownOp = (PopFuncWithUpDownOp) populationFunction;
            UpDownOperator upDownOperator1 = popFuncWithUpDownOp.getUpDownOperator1((Tree) value);
            UpDownOperator upDownOperator2 = popFuncWithUpDownOp.getUpDownOperator2((Tree) value);
            context.addExtraOperator(upDownOperator1);
            context.addExtraOperator(upDownOperator2);
        }

        if (populationFunction instanceof PopFuncWithUpOp) {
            PopFuncWithUpOp popFuncWithUpOp = (PopFuncWithUpOp) populationFunction;
            UpDownOperator upOperator = popFuncWithUpOp.getUpOperator((Tree) value);
            context.addExtraOperator(upOperator);
        }

        if (populationFunction instanceof PopFuncWithAVMNOp) {
            PopFuncWithAVMNOp popFuncWithAVMNOp = (PopFuncWithAVMNOp) populationFunction;
            AdaptableVarianceMultivariateNormalOperator avmnOp = popFuncWithAVMNOp.getAVMNOperator((Tree) value);
            context.addExtraOperator(avmnOp);
        }


        beastCoalescent.setInputValue("populationModel", populationFunction);
        beastCoalescent.initAndValidate();

        return beastCoalescent;
    }

    private beast.base.evolution.tree.coalescent.PopulationFunction convertSVSPopulationFunction(SVSPopulationFunction svsPopFunc, BEASTContext context) {
        PopulationFunction model = (PopulationFunction) svsPopFunc.getModels();
        Value<PopulationFunction> modelValue = new Value<>(model, null);

        Object beastObj = context.getBEASTObject(modelValue);
        if (!(beastObj instanceof beast.base.evolution.tree.coalescent.PopulationFunction)) {
            throw new IllegalArgumentException("Invalid type for model within SVSPopulationFunction: " + beastObj.getClass().getName());
        }

        return (beast.base.evolution.tree.coalescent.PopulationFunction) beastObj;
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
