package phylonco.lphybeast.spi;

import beast.base.evolution.datatype.DataType;
import jebl.evolution.sequences.SequenceType;
import lphy.core.model.Generator;
import lphy.core.model.Value;
import lphybeast.GeneratorToBEAST;
import lphybeast.ValueToBEAST;
import lphybeast.spi.LPhyBEASTExt;
import phylonco.beast.evolution.datatype.NucleotideDiploid16;
import phylonco.lphy.evolution.alignment.HaploidAlignment;
import phylonco.lphy.evolution.alignment.HomozygousAlignmentDistribution;
import phylonco.lphy.evolution.datatype.PhasedGenotype;
import phylonco.lphy.evolution.datatype.PhasedGenotypeFunction;
import phylonco.lphybeast.tobeast.generators.*;

import java.util.ArrayList;
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
        return new ArrayList<>();
    }

    @Override
    public List<Class<? extends GeneratorToBEAST>> getGeneratorToBEASTs() {
        return Arrays.asList( GT16ErrorModelToBEAST.class,
                GT16ToBEAST.class, GTUnphaseToBEAST.class,
                GompertzToBEAST.class, LogisticToBEAST.class, PopulationFunctionCoalescentToBEAST.class);
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
                HaploidAlignment.class);
    }

    @Override
    public List<Class<? extends Value>> getExcludedValue() {
        return new ArrayList<>();
    }


}
