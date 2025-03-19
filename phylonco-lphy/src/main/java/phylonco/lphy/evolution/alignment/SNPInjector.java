package phylonco.lphy.evolution.alignment;

import jebl.evolution.sequences.Nucleotides;
import lphy.base.evolution.SNPSampler;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.base.evolution.datatype.Variant;
import lphy.base.function.io.ReaderConst;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import phylonco.lphy.evolution.datatype.PhasedGenotype;

import java.util.ArrayList;
import java.util.List;

import static lphy.base.evolution.SNPSampler.getAmbiguousStateIndex;
import static phylonco.lphy.evolution.datatype.PhasedGenotype.getPhasedGenotypeIndex;

public class SNPInjector extends DeterministicFunction<Alignment> {
    public final String SNPName = "snp";
    public SNPInjector(@ParameterInfo(name = ReaderConst.ALIGNMENT, description = "one sequence alignment as reference") Value<Alignment> alignment,
                       @ParameterInfo(name = SNPName, description = "the SNPs inject into the reference sequence") Value<Variant[]> snps) {
        if (alignment == null) throw new IllegalArgumentException("The alignment cannot be null");
        if (alignment.value().length() > 1) throw new IllegalArgumentException("The alignment should be one sequence alignment");
        if (snps == null) throw new IllegalArgumentException("The snps cannot be null");
        setParam(ReaderConst.ALIGNMENT, alignment);
        setParam(SNPName, snps);
    }

    @GeneratorInfo(name = "injectSNP", examples = {"injectSNP.lphy"},
        description = "Add given SNPs in given alignment. If input alignment is haploid, then automatically convert non-SNP sites homozygous.")
    @Override
    public Value<Alignment> apply() {
        // get aprameters
        Alignment alignment = getAlignment().value();
        Variant[] snps = getSnps().value();

        // initialise the output alignmet
        Alignment outAlignment = new SimpleAlignment(Taxa.createTaxa(alignment.getTaxaNames()),
                alignment.nchar(), PhasedGenotype.INSTANCE);

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
        for (Integer i: snp_positions) {
            if (alignment.getSequenceTypeStr().equals(Nucleotides.NAME)){
                int inputIndex = alignment.getState(0,i);
                inputIndex = getAmbiguousStateIndex(inputIndex);
                newIndex = getPhasedGenotypeIndex(inputIndex, inputIndex);
            } else {
                newIndex = i;
            }
            outAlignment.setState(0, i, newIndex);
        }
        return null;
    }

    private int getSNPState(Variant snp) {
        int ref = snp.getRef();
        int alt = snp.getAlt();
        String genotype = snp.getGenotype();
        int state = -1;
        if (genotype == "0|0") {
            state = getPhasedGenotypeIndex(ref,ref);
        } else if (genotype == "0|1") {
            state = getPhasedGenotypeIndex(ref,alt);
        }
        return state;
    }

    public Value<Alignment> getAlignment(){
        return getParams().get(ReaderConst.ALIGNMENT);
    }

    public Value<Variant[]> getSnps(){
        return getParams().get(SNPName);
    }

}
