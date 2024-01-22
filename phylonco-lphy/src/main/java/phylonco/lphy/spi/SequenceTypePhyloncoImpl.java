package phylonco.lphy.spi;

import jebl.evolution.sequences.SequenceType;
import lphy.base.spi.SequenceTypeBaseImpl;
import lphy.base.spi.SequenceTypeExtension;
import phylonco.lphy.evolution.datatype.PhasedGenotype;
import phylonco.lphy.evolution.datatype.UnphasedGenotype;

import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/**
 * The "Container" provider class that implements SPI
 * which include a list of {@link SequenceType} to extend.
 * It requires a public no-args constructor.
 * @author Walter Xie
 */
public class SequenceTypePhyloncoImpl implements SequenceTypeExtension {

    @Override
    public Map<String, ? extends SequenceType> declareSequenceTypes() {
        Map<String, SequenceType> dataTypeMap = new ConcurrentHashMap<>();
        dataTypeMap.put(SequenceTypeBaseImpl.sanitise(PhasedGenotype.NAME), PhasedGenotype.INSTANCE);
        dataTypeMap.put(SequenceTypeBaseImpl.sanitise(UnphasedGenotype.NAME), PhasedGenotype.INSTANCE);
        return dataTypeMap;
    }

    @Override
    public Set<SequenceType> getSequenceTypes() {
        return new HashSet<>(dataTypeMap.values());
    }


    /**
     * LPhy sequence types {@link SequenceType}
     */
    private static Map<String, SequenceType> dataTypeMap;


    /**
     * Required by ServiceLoader.
     */
    public SequenceTypePhyloncoImpl() {
        dataTypeMap = new ConcurrentHashMap<>();
    }

    @Override
    public void register() {
        // sequence types
        Map<String, ? extends SequenceType> newTypes = declareSequenceTypes();
        addSequenceTypes(newTypes, dataTypeMap, "Phylonco sequence types : ");
    }

    public String getExtensionName() {
        return "Phylonco sequence types";
    }

}
