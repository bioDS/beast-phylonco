package phylonco.lphy.function;

import lphy.base.evolution.VCFLines;
import lphy.base.evolution.VCFSite;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.function.io.ReaderConst;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import phylonco.lphy.evolution.datatype.PhasedGenotype;

import java.util.ArrayList;
import java.util.List;

import static phylonco.lphy.evolution.datatype.PhasedGenotype.getNucleotideIndex;

public class writeHeterozygousSites extends DeterministicFunction<VCFLines> {
    public writeHeterozygousSites(@ParameterInfo(name = ReaderConst.ALIGNMENT, description = "the diploid alignment to identify heterozygous sites from") Value<Alignment> alignment) {
        setParam(ReaderConst.ALIGNMENT, alignment);
        if (alignment == null){
            throw new IllegalArgumentException("the alignment cannot be null!");
        }
        if (! PhasedGenotype.NAME.equals(alignment.value().getSequenceType().getName()))
            throw new IllegalArgumentException("Must be phased genotype alignment!");
    }

    @GeneratorInfo(name = "heterozygotes", examples = {"heterozygousConvert.lphy"},
            description = "write heterozygous sites in the diploid alignment in a vcf file")
    @Override
    public Value<VCFLines> apply() {
        Alignment alignment = getAlignment().value();

        // get taxa names
        List<String> taxa = new ArrayList<>();
        for (int s = 0; s < alignment.ntaxa(); s++) {
            taxa.add(alignment.getTaxonName(s));
        }
        String[] taxaNames = taxa.toArray(new String[0]);

        List<VCFSite> sites = new ArrayList<>();

        for (int i = 0; i < alignment.ntaxa(); i++) {
            for (int j = 0; j < alignment.nchar(); j++) {
                int stateIndex = alignment.getState(i,j);
                if (stateIndex == 0 || stateIndex == 5 || stateIndex == 10 || stateIndex == 15) {
                    continue;
                } else {
                    int ref = getNucleotideIndex(stateIndex)[0];
                    int alt = getNucleotideIndex(stateIndex)[1];
                    String genotype = getGenotype(ref, alt);
                    VCFSite site = new VCFSite(taxaNames[i], j, ref, alt, genotype);

                    sites.add(site);
                }
            }
        }

        VCFLines vcfLines = new VCFLines(sites);

        return new Value<>(null, vcfLines, this);
    }

    private String getGenotype(int ref, int alt) {
        String genotype = "";
        if (ref == alt){
            genotype = "0|0";
        } else if (ref !=alt){
            genotype = "0|1";
        }
        return genotype;
    }

    public Value<Alignment> getAlignment() {
        return getParams().get(ReaderConst.ALIGNMENT);
    }

}
