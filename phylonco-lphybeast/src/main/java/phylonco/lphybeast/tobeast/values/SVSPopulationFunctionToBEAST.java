package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.IntegerParameter;
import lphy.base.evolution.coalescent.PopulationFunction;
import lphy.base.evolution.coalescent.populationmodel.SVSFunction;
import lphy.base.evolution.coalescent.populationmodel.SVSPopulationFunction;
import lphy.core.model.GraphicalModelNode;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.StochasticVariableSelection;

import java.util.ArrayList;

public class SVSPopulationFunctionToBEAST implements ValueToBEAST<SVSPopulationFunction, StochasticVariableSelection> {

    @Override
    public StochasticVariableSelection valueToBEAST(Value<SVSPopulationFunction> lphyPopFuncVal, BEASTContext context) {

        StochasticVariableSelection beastPopFunc = new StochasticVariableSelection();

        SVSFunction gen = (SVSFunction) lphyPopFuncVal.getGenerator();

        IntegerParameter indicatorParam = context.getAsIntegerParameter(gen.getIndicator());


        Value<PopulationFunction[]> modelsValue = gen.getModels();
        Object[] modelsObjArray = modelsValue.value();

        PopulationFunction[] modelsArray = new PopulationFunction[modelsObjArray.length];
        for (int i = 0; i < modelsArray.length; i++) {
            modelsArray[i] = (PopulationFunction) modelsObjArray[i];
        }

        Integer lowerBound = 0;
        Integer upperBound = modelsArray.length - 1;

        indicatorParam.setInputValue("lower", lowerBound);
        indicatorParam.setInputValue("upper", upperBound);

        int currentValue = indicatorParam.getValue();
        if (currentValue < lowerBound) {
            indicatorParam.setValue(lowerBound);
        } else if (currentValue > upperBound) {
            indicatorParam.setValue(upperBound);
        }

        ArrayList<beast.base.evolution.tree.coalescent.PopulationFunction> modelFuncs = new ArrayList<>(modelsArray.length);

        for (int i = 0; i < modelsArray.length; i++) {
            String modelName = gen.MODELS_PARAM_NAME;
            GraphicalModelNode node = (GraphicalModelNode) (gen.getParams().get(modelName).getInputs().get(0));
            Object contextObj = context.getBEASTObject((Value) node.getInputs().get(i));
            beast.base.evolution.tree.coalescent.PopulationFunction beastPopFunction;
            beastPopFunction = (beast.base.evolution.tree.coalescent.PopulationFunction) contextObj;
            modelFuncs.add(i, beastPopFunction);
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
