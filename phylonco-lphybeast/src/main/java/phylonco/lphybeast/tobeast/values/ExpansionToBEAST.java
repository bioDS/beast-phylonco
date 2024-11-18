package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.ExpansionPopulation;
import lphy.base.evolution.coalescent.populationmodel.ExpansionPopulationFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.ExpansionGrowth;

public class ExpansionToBEAST implements ValueToBEAST<ExpansionPopulation, ExpansionGrowth> {

    @Override
    public ExpansionGrowth valueToBEAST(Value<ExpansionPopulation> lphyPopFuncVal, BEASTContext context) {

        // Get parameters from the LPhy value generator
        ExpansionPopulationFunction gen = (ExpansionPopulationFunction) lphyPopFuncVal.getGenerator();

        RealParameter tauParam = context.getAsRealParameter(gen.getTau());
        RealParameter rParam = context.getAsRealParameter(gen.getR());
        RealParameter ncParam = context.getAsRealParameter(gen.getNC());
        RealParameter xParam = context.getAsRealParameter(gen.getX());  // Use x instead of N0

        // Initialize ExpansionGrowth instance
        ExpansionGrowth beastPopFunc = new ExpansionGrowth();

        // Set ExpansionGrowth parameters
        beastPopFunc.setInputValue("r", rParam);
        beastPopFunc.setInputValue("NC", ncParam);
        beastPopFunc.setInputValue("tau", tauParam);
        beastPopFunc.setInputValue("x", xParam);  // Pass x as a parameter to ExpansionGrowth

        // Initialize and validate the model
        beastPopFunc.initAndValidate();

        // Set the ID to match the LPhy value
        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        return beastPopFunc;
    }

    @Override
    public Class<ExpansionPopulation> getValueClass() {
        return ExpansionPopulation.class;
    }

    @Override
    public Class<ExpansionGrowth> getBEASTClass() {
        return ExpansionGrowth.class;
    }
}
