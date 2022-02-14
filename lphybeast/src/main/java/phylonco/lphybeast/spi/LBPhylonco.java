package phylonco.lphybeast.spi;

import beast.evolution.datatype.DataType;
import jebl.evolution.sequences.SequenceType;
import lphy.graphicalModel.Generator;
import lphy.graphicalModel.Value;
import lphybeast.GeneratorToBEAST;
import lphybeast.ValueToBEAST;
import lphybeast.spi.LPhyBEASTExt;
import phylonco.lphy.evolution.datatype.PhasedGenotype;
import phylonco.lphy.evolution.datatype.PhasedGenotypeFunction;
import phylonco.lphybeast.tobeast.generators.GT16ErrorModelToBEAST;
import phylonco.lphybeast.tobeast.generators.GT16ToBEAST;
import phylonco.lphybeast.tobeast.generators.GTUnphaseToBEAST;

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
 * It requires a public no-args constructor.
 *
 * @author Walter Xie
 */
public class LBPhylonco implements LPhyBEASTExt {

    @Override
    public List<Class<? extends ValueToBEAST>> getValuesToBEASTs() {
        return Arrays.asList(new Class[0]);
    }

    @Override
    public List<Class<? extends GeneratorToBEAST>> getGeneratorToBEASTs() {
        return Arrays.asList( GT16ErrorModelToBEAST.class,
                GT16ToBEAST.class, GTUnphaseToBEAST.class );
    }

    @Override
    public Map<SequenceType, DataType> getDataTypeMap() {
        Map<SequenceType, DataType> dataTypeMap = new ConcurrentHashMap<>();
        dataTypeMap.put(PhasedGenotype.INSTANCE, new beast.evolution.datatype.NucleotideDiploid16());
        return dataTypeMap;
    }

    @Override
    public List<Class<? extends Generator>> getExcludedGenerator() {
        return Arrays.asList(PhasedGenotypeFunction.class);
    }

    @Override
    public List<Class<? extends Value>> getExcludedValue() {
        return new ArrayList<>();
    }


}
