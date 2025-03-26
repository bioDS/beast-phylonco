package phylonco.lphy.spi;

import lphy.base.spi.LPhyBaseImpl;
import lphy.core.model.BasicFunction;
import lphy.core.model.GenerativeDistribution;
import phylonco.lphy.evolution.alignment.*;
import phylonco.lphy.evolution.datatype.PhasedGenotypeFunction;
import phylonco.lphy.evolution.readcountmodel.CoverageModel;
import phylonco.lphy.evolution.readcountmodel.PloidyModel;
import phylonco.lphy.evolution.readcountmodel.ReadCountModel;
import phylonco.lphy.evolution.substitutionmodel.GT16;

import java.util.Arrays;
import java.util.List;

/**
 *
 * The provider of SPI which is an implementation of a service.
 * It requires a public no-args constructor.
 * @author Walter Xie
 */
public class PhyloncoImpl extends LPhyBaseImpl {

    /**
     * Required by ServiceLoader.
     */
    public PhyloncoImpl() {
        //TODO print package or classes info here?
    }

    @Override
    public List<Class<? extends GenerativeDistribution>> declareDistributions() {
        return Arrays.asList(
                GT16ErrorModel.class,
                HomozygousAlignmentDistribution.class,
                HeterozygousMutateAlignment.class,
                ReadCountModel.class,
                PloidyModel.class,
                CoverageModel.class);
    }

    @Override
    public List<Class<? extends BasicFunction>> declareFunctions() {
        return Arrays.asList(
                GT16.class,
                PhasedGenotypeFunction.class, UnphaseGenotypeAlignment.class,
                HaploidAlignment.class, SNPInjector.class
        );
    }

    public String getExtensionName() {
        return "Phylonco lphy library";
    }
}
