package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.ExpansionPopulation;
import lphy.base.evolution.coalescent.populationmodel.ExpansionPopulationFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.ExpansionGrowth;

/**
 * ExpansionToBEAST converts an LPhy ExpansionPopulation model (with optional I_na logic)
 * to a BEAST ExpansionGrowth model that also supports an optional NA usage.
 */
public class ExpansionToBEAST implements ValueToBEAST<ExpansionPopulation, ExpansionGrowth> {

    // If you want to fix the I_na bounds:
    private static final int I_NA_LOWER_BOUND = 0;
    private static final int I_NA_UPPER_BOUND = 1;

    @Override
    public ExpansionGrowth valueToBEAST(Value<ExpansionPopulation> lphyPopFuncVal, BEASTContext context) {

        ExpansionPopulationFunction gen = (ExpansionPopulationFunction) lphyPopFuncVal.getGenerator();

        RealParameter NCParam = context.getAsRealParameter(gen.getNC());
        RealParameter NAParam = context.getAsRealParameter(gen.getNA());
        RealParameter rParam  = context.getAsRealParameter(gen.getR());
        RealParameter xParam  = context.getAsRealParameter(gen.getX());

        IntegerParameter iNaParam = null;
        if (gen.getI_na() != null && gen.getI_na().value() != null) {
            iNaParam = context.getAsIntegerParameter(gen.getI_na());
            setIntegerParamBounds(iNaParam, I_NA_LOWER_BOUND, I_NA_UPPER_BOUND);
        }

        ExpansionGrowth beastPopFunc = new ExpansionGrowth();

        beastPopFunc.setInputValue("NC", NCParam);
        beastPopFunc.setInputValue("NA", NAParam);
        beastPopFunc.setInputValue("r",  rParam);
        beastPopFunc.setInputValue("x",  xParam);

        if (iNaParam != null) {
            beastPopFunc.setInputValue("I_na", iNaParam);
        }


        beastPopFunc.initAndValidate();

        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        return beastPopFunc;
    }

    /**
     * Helper method to ensure iNaParam is between 0 and 1
     */
    private void setIntegerParamBounds(IntegerParameter param, int lowerBound, int upperBound) {
        param.setInputValue("lower", lowerBound);
        param.setInputValue("upper", upperBound);

        int currentValue = param.getValue();
        if (currentValue < lowerBound) {
            param.setValue(lowerBound);
        } else if (currentValue > upperBound) {
            param.setValue(upperBound);
        }
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
