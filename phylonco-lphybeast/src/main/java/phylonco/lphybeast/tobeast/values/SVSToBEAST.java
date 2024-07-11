package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.PopulationFunction;
import lphy.base.evolution.coalescent.populationmodel.SVS;
import lphy.base.evolution.coalescent.populationmodel.SVSFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.StochasticVariableSelection;

public class SVSToBEAST implements ValueToBEAST<SVS, StochasticVariableSelection> {

    @Override
    public StochasticVariableSelection valueToBEAST(Value<SVS> lphyPopFuncVal, BEASTContext context) {

        StochasticVariableSelection SVSPopFunc;

        SVSFunction gen = (SVSFunction) lphyPopFuncVal.getGenerator();

        // Get the indicator and models
        RealParameter indicatorParam = context.getAsRealParameter(gen.getIndicator());
        Value<PopulationFunction[]> modelsValue = gen.getModels();
        PopulationFunction[] modelsArray = modelsValue.value();

        // Store the converted models
        beast.base.evolution.tree.coalescent.PopulationFunction[] modelFuncs = new beast.base.evolution.tree.coalescent.PopulationFunction[modelsArray.length];

        for (int i =0; i < modelsArray.length; i++) {
            modelFuncs[i] = (beast.base.evolution.tree.coalescent.PopulationFunction)context.getBEASTObject((Value<PopulationFunction>)modelsArray[i]);
        }

        // Create and return the StochasticVariableSelection instance with the converted models
        SVSPopFunc = new StochasticVariableSelection();
        SVSPopFunc.setInputValue("indicator", indicatorParam);
        SVSPopFunc.setInputValue("models", modelFuncs);
        SVSPopFunc.initAndValidate();

        ValueToParameter.setID(SVSPopFunc, lphyPopFuncVal);

        return SVSPopFunc;
    }



    @Override
    public Class getValueClass() {
        return SVS.class;
    }

    @Override
    public Class<StochasticVariableSelection> getBEASTClass() {
        return StochasticVariableSelection.class;
    }
}
