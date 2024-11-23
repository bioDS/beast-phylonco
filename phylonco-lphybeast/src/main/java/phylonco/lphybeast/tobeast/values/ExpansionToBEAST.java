package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.ExpansionPopulation;
import lphy.base.evolution.coalescent.populationmodel.ExpansionPopulationFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.ExpansionGrowth;

/**
 * ExpansionPopulationToBEAST converts an LPhy ExpansionPopulation model to a BEAST ExpansionGrowth model.
 *
 * <p>This class implements the ValueToBEAST interface to facilitate the translation of population
 * model parameters from the LPhy framework to the BEAST framework.</p>
 */
public class ExpansionToBEAST implements ValueToBEAST<ExpansionPopulation, ExpansionGrowth> {

    /**
     * Converts an LPhy ExpansionPopulation value to a BEAST ExpansionGrowth instance.
     *
     * @param lphyPopFuncVal The LPhy ExpansionPopulation value.
     * @param context        The BEAST context for parameter retrieval.
     * @return A BEAST ExpansionGrowth instance representing the same population model.
     */
    @Override
    public ExpansionGrowth valueToBEAST(Value<ExpansionPopulation> lphyPopFuncVal, BEASTContext context) {

        // Retrieve the ExpansionPopulationFunction generator from the LPhy value
        ExpansionPopulationFunction gen = (ExpansionPopulationFunction) lphyPopFuncVal.getGenerator();

        // Extract parameters from the generator using the BEAST context
        RealParameter NCParam = context.getAsRealParameter(gen.getNC());
        RealParameter NAParam = context.getAsRealParameter(gen.getNA());
        RealParameter rParam = context.getAsRealParameter(gen.getR());
        RealParameter xParam = context.getAsRealParameter(gen.getX());

        // Initialize the BEAST ExpansionGrowth instance
        ExpansionGrowth beastPopFunc = new ExpansionGrowth();

        // Set the ExpansionGrowth parameters with the corresponding BEAST RealParameters
        beastPopFunc.setInputValue("NC", NCParam);
        beastPopFunc.setInputValue("NA", NAParam);
        beastPopFunc.setInputValue("r", rParam);
        beastPopFunc.setInputValue("x", xParam);

        // Initialize and validate the BEAST ExpansionGrowth model
        beastPopFunc.initAndValidate();

        // Assign a unique ID to the BEAST model to match the LPhy value
        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        // Return the configured BEAST ExpansionGrowth instance
        return beastPopFunc;
    }

    /**
     * Specifies the LPhy value class that this converter handles.
     *
     * @return The ExpansionPopulation class.
     */
    @Override
    public Class<ExpansionPopulation> getValueClass() {
        return ExpansionPopulation.class;
    }

    /**
     * Specifies the BEAST class that this converter produces.
     *
     * @return The ExpansionGrowth class.
     */
    @Override
    public Class<ExpansionGrowth> getBEASTClass() {
        return ExpansionGrowth.class;
    }
}
