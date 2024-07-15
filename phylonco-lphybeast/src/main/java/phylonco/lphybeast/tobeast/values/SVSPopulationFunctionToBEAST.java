package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.PopulationFunction;
import lphy.base.evolution.coalescent.populationmodel.SVSPopulationFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.StochasticVariableSelection;

public class SVSPopulationFunctionToBEAST implements ValueToBEAST<SVSPopulationFunction, StochasticVariableSelection> {

    @Override
    public StochasticVariableSelection valueToBEAST(Value<SVSPopulationFunction> lphyPopFuncVal, BEASTContext context) {

        StochasticVariableSelection beastPopFunc = new StochasticVariableSelection();

        SVSPopulationFunction gen = (SVSPopulationFunction) lphyPopFuncVal.getGenerator();


        RealParameter indicatorParam = context.getAsRealParameter(gen.getIndicator());
        Value<PopulationFunction[]> modelsValue = gen.getModels();
        PopulationFunction[] modelsArray = modelsValue.value();


        beast.base.evolution.tree.coalescent.PopulationFunction[] modelFuncs = new beast.base.evolution.tree.coalescent.PopulationFunction[modelsArray.length];

        for (int i = 0; i < modelsArray.length; i++) {
            modelFuncs[i] = (beast.base.evolution.tree.coalescent.PopulationFunction) context.getBEASTObject(new Value<>(String.valueOf(modelsArray[i]), gen));
        }


        beastPopFunc.setInputValue("indicator", indicatorParam);
        beastPopFunc.setInputValue("models", modelFuncs);
        beastPopFunc.initAndValidate();

        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        return beastPopFunc;
    }

    @Override
    public Class getValueClass() {
        return SVSPopulationFunction.class;
    }

    @Override
    public Class<StochasticVariableSelection> getBEASTClass() {
        return StochasticVariableSelection.class;
    }
}
