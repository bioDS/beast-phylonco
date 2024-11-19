package phylonco.lphybeast.tobeast.values;

import beast.base.inference.parameter.RealParameter;
import lphy.base.evolution.coalescent.populationmodel.ExponentialPopulation;
import lphy.base.evolution.coalescent.populationmodel.ExponentialPopulationFunction;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.ValueToBEAST;
import lphybeast.tobeast.values.ValueToParameter;
import phylonco.beast.evolution.populationmodel.ExponentialGrowth;

public class ExponentialToBEAST implements ValueToBEAST<ExponentialPopulation, ExponentialGrowth> {

    @Override
    public ExponentialGrowth valueToBEAST(Value<ExponentialPopulation> lphyPopFuncVal, BEASTContext context) {

        ExponentialGrowth beastPopFunc;

        // 获取 ExponentialPopulationFunction 的实例
        ExponentialPopulationFunction gen = (ExponentialPopulationFunction) lphyPopFuncVal.getGenerator();

        // 获取 GrowthRate 和 N0 参数
        RealParameter GrowthRateParam = context.getAsRealParameter(gen.getGrowthRate());
        RealParameter N0Param = context.getAsRealParameter(gen.getN0());

        beastPopFunc = new ExponentialGrowth();

        // 设置 GrowthRate 和 N0 输入
        beastPopFunc.setInputValue("GrowthRate", GrowthRateParam);
        beastPopFunc.setInputValue("N0", N0Param);

        // 检查是否包含 NA 参数
        Value<Double> naValue = gen.getNA();
        if (naValue != null && naValue.value() != null && naValue.value() > 0.0) {
            // 获取 NA 参数
            RealParameter NAParam = context.getAsRealParameter(gen.getNA());
            // 设置 NA 输入
            beastPopFunc.setInputValue("NA", NAParam);
        }

        // 初始化并验证 BEAST 模型
        beastPopFunc.initAndValidate();

        // 设置参数 ID
        ValueToParameter.setID(beastPopFunc, lphyPopFuncVal);

        return beastPopFunc;
    }

    @Override
    public Class getValueClass() {
        return ExponentialPopulation.class;
    }

    @Override
    public Class<ExponentialGrowth> getBEASTClass() {
        return ExponentialGrowth.class;
    }

}
