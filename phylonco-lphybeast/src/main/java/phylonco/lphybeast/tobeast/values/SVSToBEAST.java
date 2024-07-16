package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.PopulationFunction;
import lphy.base.evolution.coalescent.populationmodel.SVSFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import phylonco.beast.evolution.populationmodel.StochasticVariableSelection;

public class SVSToBEAST implements ValueToBEAST<SVSFunction, StochasticVariableSelection> {

    @Override
    public StochasticVariableSelection valueToBEAST(Value<SVSFunction> lphyPopFuncVal, BEASTContext context) {

        SVSFunction gen = (SVSFunction) lphyPopFuncVal.getGenerator();

        // Get the indicator and models
        RealParameter indicatorParam = context.getAsRealParameter(gen.getIndicator());
        PopulationFunction[] modelsArray = gen.getModels().value();

        // Store the converted models
        beast.base.evolution.tree.coalescent.PopulationFunction[] modelFuncs = new beast.base.evolution.tree.coalescent.PopulationFunction[modelsArray.length];

        for (int i = 0; i < modelsArray.length; i++) {
            modelFuncs[i] = (beast.base.evolution.tree.coalescent.PopulationFunction) context.getBEASTObject(new Value<>(null, modelsArray[i]));
        }

        // Create and return the StochasticVariableSelection instance with the converted models
        StochasticVariableSelection stochasticGrowth = new StochasticVariableSelection();
        stochasticGrowth.setInputValue("indicator", indicatorParam);
        stochasticGrowth.setInputValue("models", modelFuncs);
        stochasticGrowth.initAndValidate();

        return stochasticGrowth;
    }

    @Override
    public Class<SVSFunction> getValueClass() {
        return SVSFunction.class;
    }

    @Override
    public Class<StochasticVariableSelection> getBEASTClass() {
        return StochasticVariableSelection.class;
    }
}