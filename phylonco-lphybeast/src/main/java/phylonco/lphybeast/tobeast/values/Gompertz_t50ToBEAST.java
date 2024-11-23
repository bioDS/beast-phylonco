package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulationFunction_t50;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulation_t50;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.GompertzGrowth_t50;

public class Gompertz_t50ToBEAST implements ValueToBEAST<GompertzPopulation_t50, GompertzGrowth_t50> {

    @Override
    public GompertzGrowth_t50 valueToBEAST(Value<GompertzPopulation_t50> lphyPopFuncVal, BEASTContext context) {

        GompertzGrowth_t50 beastPopFunc = new GompertzGrowth_t50();

        GompertzPopulationFunction_t50 gen = (GompertzPopulationFunction_t50) lphyPopFuncVal.getGenerator();

        RealParameter t50Param = context.getAsRealParameter(gen.getT50());
        RealParameter bParam = context.getAsRealParameter(gen.getB());
        RealParameter NInfinityParam = context.getAsRealParameter(gen.getNInfinity());

        RealParameter NAParam = null;
        if (gen.getNA() != null && gen.getNA().value() != null) {
            NAParam = context.getAsRealParameter(gen.getNA());
        }

        beastPopFunc.setInputValue("t50", t50Param);
        beastPopFunc.setInputValue("b", bParam);
        beastPopFunc.setInputValue("NInfinity", NInfinityParam);

        if (NAParam != null) {
            beastPopFunc.setInputValue("NA", NAParam);
        }

        beastPopFunc.initAndValidate();

        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        return beastPopFunc;
    }

    @Override
    public Class<GompertzPopulation_t50> getValueClass() {
        return GompertzPopulation_t50.class;
    }

    @Override
    public Class<GompertzGrowth_t50> getBEASTClass() {
        return GompertzGrowth_t50.class;
    }

}
