package phylonco.lphy.evolution.alignment;

import jebl.evolution.sequences.Nucleotides;
import jebl.evolution.sequences.SequenceType;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.base.evolution.datatype.DataType;
import lphy.base.evolution.datatype.Variant;
import lphy.base.evolution.likelihood.AbstractPhyloCTMC;
import lphy.base.function.io.ReaderConst;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import phylonco.lphy.evolution.datatype.PhasedGenotype;
import phylonco.lphy.evolution.datatype.UnphasedGenotype;

import java.util.ArrayList;
import java.util.List;

import static phylonco.lphy.evolution.alignment.HomozygousAlignmentDistribution.getHomozygousState;
import static phylonco.lphy.evolution.datatype.PhasedGenotype.getPhasedGenotypeIndex;
import static phylonco.lphy.evolution.datatype.UnphasedGenotype.getUnphasedGenotypeIndex;

public class SNPInjector extends DeterministicFunction<Alignment> {
    public final String SNPName = "snp";
    protected Value<DataType> dataType;
    public SNPInjector(@ParameterInfo(name = ReaderConst.ALIGNMENT, description = "one sequence alignment as reference") Value<Alignment> alignment,
                       @ParameterInfo(name = AbstractPhyloCTMC.dataTypeParamName, description = "the data type used for simulations, default to nucleotide",
                               narrativeName = "data type used for simulations", optional = true) Value<DataType> dataType,
                       @ParameterInfo(name = SNPName, description = "the SNPs inject into the reference sequence") Value<Variant[]> snps) {
        if (alignment == null) throw new IllegalArgumentException("The alignment cannot be null");
        if (alignment.value().length() > 1) throw new IllegalArgumentException("The alignment should be one sequence alignment");
        if (snps == null) throw new IllegalArgumentException("The snps cannot be null");
        if (dataType == null) throw new IllegalArgumentException("The data type cannot be null");
        setParam(ReaderConst.ALIGNMENT, alignment);
        setParam(SNPName, snps);
        this.dataType = dataType;
    }

    @GeneratorInfo(name = "injectSNP", examples = {"SNPInjection.lphy"},
        description = "Add given SNPs in given alignment. If input alignment is haploid, then automatically convert non-SNP sites homozygous.")
    @Override
    public Value<Alignment> apply() {
        // get aprameters
        Alignment alignment = getAlignment().value();
        Variant[] snps = getSNPs().value();

        // initialise the output alignment
        Alignment outAlignment;
        if (dataType.value() instanceof PhasedGenotype) {
            outAlignment = new SimpleAlignment(Taxa.createTaxa(alignment.getTaxaNames()),
                    alignment.nchar(), PhasedGenotype.INSTANCE);
        } else if (dataType.value() instanceof UnphasedGenotype) {
            outAlignment = new SimpleAlignment(Taxa.createTaxa(alignment.getTaxaNames()),
                    alignment.nchar(), UnphasedGenotype.INSTANCE);
        } else {
            throw new IllegalArgumentException(
                    "Unsupported data type: " + dataType.value().getClass().getName());
        }


        // map the alignment
        List<Integer> snp_positions = new ArrayList<>();
        for (Variant snp : snps) {
            int position = snp.getPosition();
            snp_positions.add(position);
            int state = getSNPState(snp);
            outAlignment.setState(0, position, state);
        }

        // fill other positions
        int newIndex = -1;
        for (int i = 0; i < alignment.nchar(); i++){
            if (! snp_positions.contains(i)) {
                if (alignment.getSequenceTypeStr().equals(Nucleotides.NAME)){
                    newIndex = getHomozygousState(alignment, 0 ,i, dataType.value());
                } else {
                    newIndex = i;
                }
                outAlignment.setState(0, i, newIndex);
            }
        }

        return new Value<>("A", outAlignment, this);
    }

    private int getSNPState(Variant snp) {
        int ref = snp.getRef();
        int alt = snp.getAlt();
        String genotype = snp.getGenotype();
        int state = -1;
        if (dataType.value() instanceof PhasedGenotype) {
            if (genotype == "0|0") {
                state = getPhasedGenotypeIndex(ref,ref);
            } else if (genotype == "0|1") {
                state = getPhasedGenotypeIndex(ref,alt);
            } else if (genotype == "1|1") {
                state = getPhasedGenotypeIndex(alt,alt);
            }
        } else if (dataType.value() instanceof UnphasedGenotype) {
            if (genotype == "0|0") {
                state = getUnphasedGenotypeIndex(ref,ref);
            } else if (genotype == "0|1") {
                state = getUnphasedGenotypeIndex(ref,alt);
            } else if (genotype == "1|1") {
                state = getUnphasedGenotypeIndex(alt,alt);
            }
        } else {
            throw new IllegalArgumentException(
                    "Unsupported data type: " + dataType.value().getClass().getName());
        }

        return state;
    }

    public Value<Alignment> getAlignment(){
        return getParams().get(ReaderConst.ALIGNMENT);
    }

    public Value<Variant[]> getSNPs(){
        return getParams().get(SNPName);
    }

}
