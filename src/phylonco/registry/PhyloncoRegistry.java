package phylonco.registry;

import beast.evolution.datatype.DataType;
import jebl.evolution.sequences.SequenceType;
import lphy.core.LPhyParser;
import lphy.evolution.datatype.PhasedGenotype;
import lphybeast.registry.ClassesRegistry;
import phylonco.tobeast.generators.GT16ErrorModelToBEAST;
import phylonco.tobeast.generators.GT16ToBEAST;
import phylonco.tobeast.generators.GTUnphaseToBEAST;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Registry class name needs to be ended with the keyword "Registry" (case insensitive),
 * in order to be registered by {@link lphybeast.BEASTContext#BEASTContext(LPhyParser)}.
 *
 * @author Walter Xie
 */
public class PhyloncoRegistry implements ClassesRegistry {

    private final Class<?>[] generatorToBEASTs = {
            GT16ErrorModelToBEAST.class,
            GT16ToBEAST.class,
            GTUnphaseToBEAST.class
    };

    @Override
    public Class<?>[] getValuesToBEASTs() {
        return new Class[0];
    }

    @Override
    public Class<?>[] getGeneratorToBEASTs() {
        return generatorToBEASTs;
    }

    @Override
    public Map<SequenceType, DataType> getDataTypeMap() {
        Map<SequenceType, DataType> dataTypeMap = new ConcurrentHashMap<>();
        dataTypeMap.put(PhasedGenotype.INSTANCE, new beast.evolution.datatype.NucleotideDiploid16());
        return dataTypeMap;
    }
}
