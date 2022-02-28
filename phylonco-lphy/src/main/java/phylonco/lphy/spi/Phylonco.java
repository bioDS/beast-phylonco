package phylonco.lphy.spi;

import jebl.evolution.sequences.SequenceType;
import lphy.evolution.datatype.SequenceTypeFactory;
import lphy.graphicalModel.Func;
import lphy.graphicalModel.GenerativeDistribution;
import lphy.spi.LPhyExtension;
import phylonco.lphy.evolution.alignment.GT16ErrorModel;
import phylonco.lphy.evolution.alignment.UnphaseGenotypeAlignment;
import phylonco.lphy.evolution.datatype.PhasedGenotype;
import phylonco.lphy.evolution.datatype.PhasedGenotypeFunction;
import phylonco.lphy.evolution.datatype.UnphasedGenotype;
import phylonco.lphy.evolution.substitutionmodel.GT16;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

/**
 *
 * The provider of SPI which is an implementation of a service.
 * It requires a public no-args constructor.
 * @author Walter Xie
 */
public class Phylonco implements LPhyExtension {

    /**
     * Required by ServiceLoader.
     */
    public Phylonco() {
        //TODO print package or classes info here?
    }

    @Override
    public List<Class<? extends GenerativeDistribution>> getDistributions() {
        return Arrays.asList( GT16ErrorModel.class );
    }

    @Override
    public List<Class<? extends Func>> getFunctions() {
        return Arrays.asList( GT16.class, PhasedGenotypeFunction.class, UnphaseGenotypeAlignment.class);
    }

    @Override
    public Map<String, ? extends SequenceType> getSequenceTypes() {
        Map<String, SequenceType> dataTypeMap = new ConcurrentHashMap<>();
        dataTypeMap.put(SequenceTypeFactory.sanitise(PhasedGenotype.NAME), PhasedGenotype.INSTANCE);
        dataTypeMap.put(SequenceTypeFactory.sanitise(UnphasedGenotype.NAME), PhasedGenotype.INSTANCE);
        return dataTypeMap;
    }
}
