package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.Cons_Exp_ConsPopulation;
import lphy.base.evolution.coalescent.populationmodel.Cons_Exp_ConsPopulationFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.Cons_Exp_ConsGrowth;

public class Cons_Exp_ConsToBEAST implements ValueToBEAST<Cons_Exp_ConsPopulation, Cons_Exp_ConsGrowth> {

    @Override
    public Cons_Exp_ConsGrowth valueToBEAST(Value<Cons_Exp_ConsPopulation> lphyPopFuncVal, BEASTContext context) {

        // Get parameters from the LPhy value generator
        Cons_Exp_ConsPopulationFunction gen = (Cons_Exp_ConsPopulationFunction) lphyPopFuncVal.getGenerator();

        RealParameter tauParam = context.getAsRealParameter(gen.getTau());
        RealParameter rParam = context.getAsRealParameter(gen.getR());
        RealParameter ncParam = context.getAsRealParameter(gen.getNC());
        RealParameter xParam = context.getAsRealParameter(gen.getX());  // Use x instead of N0

        // Initialize Cons_Exp_ConsGrowth instance
        Cons_Exp_ConsGrowth beastPopFunc = new Cons_Exp_ConsGrowth();

        // Set Cons_Exp_ConsGrowth parameters
        beastPopFunc.setInputValue("r", rParam);
        beastPopFunc.setInputValue("NC", ncParam);
        beastPopFunc.setInputValue("tau", tauParam);
        beastPopFunc.setInputValue("x", xParam);  // Pass x as a parameter to Cons_Exp_ConsGrowth

        // Initialize and validate the model
        beastPopFunc.initAndValidate();

        // Set the ID to match the LPhy value
        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        return beastPopFunc;
    }

    @Override
    public Class<Cons_Exp_ConsPopulation> getValueClass() {
        return Cons_Exp_ConsPopulation.class;
    }

    @Override
    public Class<Cons_Exp_ConsGrowth> getBEASTClass() {
        return Cons_Exp_ConsGrowth.class;
    }
}