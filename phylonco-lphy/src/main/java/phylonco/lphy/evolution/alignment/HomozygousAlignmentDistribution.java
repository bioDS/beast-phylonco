package phylonco.lphy.evolution.alignment;

import lphy.base.distribution.ParametricDistribution;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.base.function.io.ReaderConst;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.random.RandomGenerator;
import phylonco.lphy.evolution.datatype.PhasedGenotype;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import static lphy.base.evolution.SNPSampler.getAmbiguousStateIndex;
import static phylonco.lphy.evolution.datatype.PhasedGenotype.getPhasedGenotypeIndex;

public class HomozygousAlignmentDistribution extends ParametricDistribution<Alignment> {
    private Value<Alignment> alignmentValue;

    public HomozygousAlignmentDistribution(@ParameterInfo(name = ReaderConst.ALIGNMENT,
            description = "Convert the input haploid alignment into genotype alignment (homozygous).")
                                           Value<Alignment> alignmentValue) {
        if (alignmentValue == null) throw new IllegalArgumentException("The alignment can't be null!");
        this.alignmentValue = alignmentValue;
    }

    @Override
    protected void constructDistribution(RandomGenerator random) {
    }

    @GeneratorInfo(name = "Homozygous", description = "Convert the haploid to homozygous alignment. The transformation is deterministic" +
            " when there are no ambiguities in the input alignment. If there are ambiguous states or gaps in the sequence," +
            "states will be chosen randomly among the possible states and give a homozygous alignment.")
    @Override
    public RandomVariable<Alignment> sample() {
        // get the original seq
        Alignment originalAlignment = ((Value<Alignment>) getParams().get(ReaderConst.ALIGNMENT)).value();

        // obtain the taxa names
        List<String> taxa = new ArrayList<>();
        for (int s = 0; s < originalAlignment.ntaxa(); s++) {
            taxa.add(originalAlignment.getTaxonName(s));
        }
        String[] taxaNames = taxa.toArray(new String[0]);

        // initialise the new alignment
        Alignment genotypeAlignment = new SimpleAlignment(Taxa.createTaxa(taxaNames),
                originalAlignment.nchar(), PhasedGenotype.INSTANCE);

        // set the alignment
        for (int i = 0; i < genotypeAlignment.ntaxa(); i++) {
            for (int j = 0; j < genotypeAlignment.nchar(); j++) {
                int index = getHomozygousState(originalAlignment, i, j);

                // map the new alignment states
                genotypeAlignment.setState(i, j, index);
            }
        }
        return new RandomVariable<>(null, genotypeAlignment, this);
    }

    public static int getHomozygousState(Alignment originalAlignment, int i, int j) {
        // get the nucleotide index of each site
        int originalStateIndex = originalAlignment.getState(i, j);

        // get the certain nucleotide index for each site
        int stateIndex = getAmbiguousStateIndex(originalStateIndex);

        // convert the nucleotide states into phased genotypes
        int index = getPhasedGenotypeIndex(stateIndex, stateIndex);
        return index;
    }

    @Override
    public Map<String, Value> getParams() {
        return new TreeMap<>() {{
            put(ReaderConst.ALIGNMENT, alignmentValue);
        }};
    }

    public void setParam(String paramName, Value value) {
        if (paramName.equals(ReaderConst.ALIGNMENT)) alignmentValue = value;
    }
}
