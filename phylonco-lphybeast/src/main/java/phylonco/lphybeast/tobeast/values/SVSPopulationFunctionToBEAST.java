package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.PopulationFunction;
import lphy.base.evolution.coalescent.populationmodel.SVSFunction;
import lphy.base.evolution.coalescent.populationmodel.SVSPopulationFunction;
import lphy.core.model.GraphicalModelNode;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.StochasticVariableSelection;

public class SVSPopulationFunctionToBEAST implements ValueToBEAST<SVSPopulationFunction, StochasticVariableSelection> {

    @Override
    public StochasticVariableSelection valueToBEAST(Value<SVSPopulationFunction> lphyPopFuncVal, BEASTContext context) {

        StochasticVariableSelection beastPopFunc = new StochasticVariableSelection();

        SVSFunction gen = (SVSFunction) lphyPopFuncVal.getGenerator();


        RealParameter indicatorParam = context.getAsRealParameter(gen.getIndicator());
        Value<PopulationFunction[]> modelsValue = gen.getModels();
        Object[] modelsObjArray = modelsValue.value();

        PopulationFunction[] modelsArray = new PopulationFunction[modelsObjArray.length];
        for (int i = 0; i < modelsArray.length; i++) {
            modelsArray[i] = (PopulationFunction) modelsObjArray[i];
        }

        beast.base.evolution.tree.coalescent.PopulationFunction[] modelFuncs = new beast.base.evolution.tree.coalescent.PopulationFunction[modelsArray.length];

        for (int i = 0; i < modelsArray.length; i++) {
            String modelName = gen.MODELS_PARAM_NAME;
            GraphicalModelNode node = (GraphicalModelNode) (gen.getParams().get(modelName).getInputs().get(0));
            Object contextObj = context.getBEASTObject((Value) node.getInputs().get(i));
            beast.base.evolution.tree.coalescent.PopulationFunction beastPopFunction;
            beastPopFunction = (beast.base.evolution.tree.coalescent.PopulationFunction) contextObj;
            modelFuncs[i] = beastPopFunction;
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
