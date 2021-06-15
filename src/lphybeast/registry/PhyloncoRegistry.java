package lphybeast.registry;

import beast.evolution.datatype.DataType;
import jebl.evolution.sequences.SequenceType;
import lphy.evolution.datatype.PhasedGenotype;
import lphybeast.phylonco.generators.GT16ErrorModelToBEAST;
import lphybeast.phylonco.generators.GT16ToBEAST;
import lphybeast.phylonco.generators.GTUnphaseToBEAST;
import lphybeast.registry.ClassesRegistry;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
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
