package phylonco.lphybeast.spi;

import beast.base.evolution.datatype.DataType;
import jebl.evolution.sequences.SequenceType;
import lphy.base.distribution.UniformDiscrete;

import lphy.base.evolution.tree.*;
import lphy.base.function.Difference;
import lphy.base.function.Union;
import lphy.base.function.io.ReadTrees;
import lphy.core.model.Generator;
import lphybeast.GeneratorToBEAST;
import lphybeast.ValueToBEAST;
import lphybeast.spi.LPhyBEASTExt;
import phylonco.beast.evolution.datatype.NucleotideDiploid10;
import phylonco.beast.evolution.datatype.NucleotideDiploid16;
import phylonco.lphy.evolution.alignment.HaploidAlignment;
import phylonco.lphy.evolution.alignment.HomozygousAlignmentDistribution;
import phylonco.lphy.evolution.alignment.SNPInjector;
import phylonco.lphy.evolution.copynumbermodel.CopyNumberBD;
import phylonco.lphy.evolution.copynumbermodel.ReadCopyProfile;
import phylonco.lphy.evolution.datatype.PhasedGenotype;
import phylonco.lphy.evolution.datatype.PhasedGenotypeFunction;
import phylonco.lphy.evolution.datatype.UnphasedGenotype;
import phylonco.lphy.evolution.datatype.UnphasedGenotypeFunction;
import phylonco.lphy.evolution.readcountmodel.CoverageModel;
import phylonco.lphy.evolution.readcountmodel.Integer2DMatrix;
import phylonco.lphy.evolution.readcountmodel.PloidyModel;
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
 *
 * @author Walter Xie
 */
public class LBPhylonco implements LPhyBEASTExt {

    @Override
    public List<Class<? extends ValueToBEAST>> getValuesToBEASTs() {
        return Arrays.asList(

                ReadCountToBEAST.class,
                //copy number model
                IntegerCharacterMatrixToBEAST.class,
                CopyNumberBDToBEAST.class
        );
    }

    @Override
    public List<Class<? extends GeneratorToBEAST>> getGeneratorToBEASTs() {
        return Arrays.asList(GT16ErrorModelToBEAST.class,
                GT16ToBEAST.class, GT10ToBEAST.class, GTUnphaseToBEAST.class,
                ReadCountModelToBEAST.class,
//                LocalClockToBeast.class//
                // copy number model
                PhyloDiscreteToBEAST.class
        );
    }

    @Override
    public Map<SequenceType, DataType> getDataTypeMap() {
        Map<SequenceType, DataType> dataTypeMap = new ConcurrentHashMap<>();
        dataTypeMap.put(PhasedGenotype.INSTANCE, new NucleotideDiploid16());
        dataTypeMap.put(UnphasedGenotype.INSTANCE, new NucleotideDiploid10());
        return dataTypeMap;
    }

    @Override
    public List<Class<? extends Generator>> getExcludedGenerator() {
        return Arrays.asList(PhasedGenotypeFunction.class, UnphasedGenotypeFunction.class, HomozygousAlignmentDistribution.class,
                HaploidAlignment.class, Difference.class, Union.class, ReadTrees.class,
                SampleBranch.class, SubstituteClade.class, SubsampledTree.class, LabelClade.class,
                MRCA.class, SNPInjector.class,
                PloidyModel.class, CoverageModel.class,
                UniformDiscrete.class,
                CopyNumberBD.class,
                ReadCopyProfile.class
        );
    }

    @Override
    public List<Class> getExcludedValueType() {
        return Arrays.asList(TimeTreeNode.class,
                Integer2DMatrix.class
        );
    }

}
