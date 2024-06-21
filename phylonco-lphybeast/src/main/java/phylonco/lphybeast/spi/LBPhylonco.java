package phylonco.lphybeast.spi;

import beast.base.evolution.datatype.DataType;
import jebl.evolution.sequences.SequenceType;
import lphy.base.evolution.branchrate.LocalClock;
import lphy.base.evolution.coalescent.PopulationFunctionCoalescent;
import lphy.base.evolution.coalescent.populationmodel.ExponentialPopulationFunction;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulationFunction_f0;
import lphy.base.evolution.coalescent.populationmodel.GompertzPopulationFunction_t50;
import lphy.base.evolution.coalescent.populationmodel.LogisticPopulationFunction;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Generator;
import lphybeast.GeneratorToBEAST;
import lphybeast.ValueToBEAST;
import lphybeast.spi.LPhyBEASTExt;
import phylonco.beast.evolution.datatype.NucleotideDiploid16;
import phylonco.lphy.evolution.alignment.HaploidAlignment;
import phylonco.lphy.evolution.alignment.HomozygousAlignmentDistribution;
import phylonco.lphy.evolution.datatype.PhasedGenotype;
import phylonco.lphy.evolution.datatype.PhasedGenotypeFunction;
import phylonco.lphybeast.tobeast.generators.GT16ErrorModelToBEAST;
import phylonco.lphybeast.tobeast.generators.GT16ToBEAST;
import phylonco.lphybeast.tobeast.generators.GTUnphaseToBEAST;
import phylonco.lphybeast.tobeast.generators.PopFuncCoalescentToBEAST;
import phylonco.lphybeast.tobeast.values.ExponentialToBEAST;
import phylonco.lphybeast.tobeast.values.Gompertz_f0ToBEAST;
import phylonco.lphybeast.tobeast.values.Gompertz_t50ToBEAST;
import phylonco.lphybeast.tobeast.values.LogisticToBEAST;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * The "Container" provider class of SPI
 * which include a list of {@link lphybeast.ValueToBEAST},
 * {@link lphybeast.GeneratorToBEAST}, and {@link DataType}
 * to extend.
 * @author Walter Xie
 */
public class LBPhylonco implements LPhyBEASTExt {

    @Override
    public List<Class<? extends ValueToBEAST>> getValuesToBEASTs() {
        return Arrays.asList(
                //PopulationFunctionToBEAST.class // TODO
                Gompertz_f0ToBEAST.class , ExponentialToBEAST.class, LogisticToBEAST.class, Gompertz_t50ToBEAST.class
        );
    }

    @Override
    public List<Class<? extends GeneratorToBEAST>> getGeneratorToBEASTs() {
        return Arrays.asList( GT16ErrorModelToBEAST.class,
                GT16ToBEAST.class, GTUnphaseToBEAST.class,
                PopFuncCoalescentToBEAST.class//, GompertzToBEAST.class
//                , LogisticToBEAST.class, PopulationFunctionCoalescentToBEAST.class
        );
    }

    @Override
    public Map<SequenceType, DataType> getDataTypeMap() {
        Map<SequenceType, DataType> dataTypeMap = new ConcurrentHashMap<>();
        dataTypeMap.put(PhasedGenotype.INSTANCE, new NucleotideDiploid16());
        return dataTypeMap;
    }

    @Override
    public List<Class<? extends Generator>> getExcludedGenerator() {
        return Arrays.asList(PhasedGenotypeFunction.class, HomozygousAlignmentDistribution.class,
                HaploidAlignment.class,
                LocalClock.class, //TODO
                PopulationFunctionCoalescent.class,
                GompertzPopulationFunction_f0.class, GompertzPopulationFunction_t50.class, ExponentialPopulationFunction.class, LogisticPopulationFunction.class);
    }

    @Override
    public List<Class> getExcludedValueType() {
        return Arrays.asList(TimeTreeNode.class);
    }

}
