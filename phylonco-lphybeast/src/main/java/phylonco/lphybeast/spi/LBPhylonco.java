package phylonco.lphybeast.spi;

import beast.base.evolution.datatype.DataType;
import jebl.evolution.sequences.SequenceType;
import lphy.base.evolution.coalescent.PopulationFunctionCoalescent;
import lphy.base.evolution.coalescent.populationmodel.*;
import lphy.base.evolution.tree.*;
import lphy.base.function.Difference;
import lphy.base.function.Union;
import lphy.base.function.io.ReadTrees;
import lphy.core.model.Generator;
import lphybeast.GeneratorToBEAST;
import lphybeast.ValueToBEAST;
import lphybeast.spi.LPhyBEASTExt;
import lphybeast.tobeast.generators.CalibratedYuleToBeast;
import phylonco.beast.evolution.datatype.NucleotideDiploid16;
import phylonco.lphy.evolution.alignment.HaploidAlignment;
import phylonco.lphy.evolution.alignment.HomozygousAlignmentDistribution;
import phylonco.lphy.evolution.datatype.PhasedGenotype;
import phylonco.lphy.evolution.datatype.PhasedGenotypeFunction;
import phylonco.lphybeast.tobeast.generators.*;
import phylonco.lphybeast.tobeast.values.*;

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
                Gompertz_f0ToBEAST.class , ExponentialToBEAST.class, LogisticToBEAST.class, Gompertz_t50ToBEAST.class, ConstantToBEAST.class
        );
    }

    @Override
    public List<Class<? extends GeneratorToBEAST>> getGeneratorToBEASTs() {
        return Arrays.asList( GT16ErrorModelToBEAST.class,
                GT16ToBEAST.class, GTUnphaseToBEAST.class,
                PopFuncCoalescentToBEAST.class,
                LocalClockToBeast.class, CalibratedYuleToBeast.class//, GompertzToBEAST.class
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
                HaploidAlignment.class, Difference.class, Union.class, ReadTrees.class,
                SampleBranch.class, SubstituteClade.class, SubsampledTree.class, LabelClade.class,
                PopulationFunctionCoalescent.class,
                ConstantPopulationFunction.class, GompertzPopulationFunction_f0.class, GompertzPopulationFunction_t50.class, ExponentialPopulationFunction.class, LogisticPopulationFunction.class);
    }

    @Override
    public List<Class> getExcludedValueType() {
        return Arrays.asList(TimeTreeNode.class);
    }

}
