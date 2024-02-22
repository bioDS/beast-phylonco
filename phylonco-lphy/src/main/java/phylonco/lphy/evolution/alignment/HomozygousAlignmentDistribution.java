package phylonco.lphy.evolution.alignment;

import lphy.base.distribution.ParametricDistribution;
import lphy.base.distribution.UniformDiscrete;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.AbstractAlignment;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.base.function.io.ReaderConst;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.random.RandomGenerator;
import phylonco.lphy.evolution.datatype.PhasedGenotype;


import java.util.*;

import static phylonco.lphy.evolution.datatype.NucleotideGenotypeHelper.getPhasedGenotypeIndex;

public class HomozygousAlignmentDistribution extends ParametricDistribution<Alignment> {
    private Value<Alignment> alignmentValue;
    public HomozygousAlignmentDistribution(@ParameterInfo(name = ReaderConst.ALIGNMENT,
            description = "Convert the input haploid alignment into genotype alignment (homozygous)." )
                                           Value<Alignment> alignmentValue){
        if (alignmentValue == null) throw new IllegalArgumentException("The alignment can't be null!");
        this.alignmentValue = alignmentValue;
    }

    @Override
    protected void constructDistribution(RandomGenerator random) {}

    @GeneratorInfo(name = "Homozygous", description = "Convert the haploid to homozygous alignment. The transformation is deterministic" +
            " when there are no ambiguities in the input alignment. If there are ambiguous states or gaps in the sequence," +
            "states will be chosen randomly among the possible states and give a homozygous alignment.")
    @Override
    public RandomVariable<Alignment> sample() {
        // get the original seq
        Alignment originalAlignment = ((Value<Alignment>) getParams().get(ReaderConst.ALIGNMENT)).value();

        // obtain the taxa names
        List<String> taxa = new ArrayList<>();
        for (int s = 0; s < originalAlignment.ntaxa(); s++){
            taxa.add(originalAlignment.getTaxonName(s));
        }
        String[] taxaNames = taxa.toArray(new String[0]);

        // initialise the new alignment
        Alignment genotypeAlignment = new SimpleAlignment(Taxa.createTaxa(taxaNames),
                originalAlignment.nchar(), PhasedGenotype.INSTANCE);

        // set the alignment
        for (int i = 0; i < genotypeAlignment.ntaxa(); i++) {
            for (int j = 0; j < genotypeAlignment.nchar(); j++) {
                // get the nucleotide index of each site
                int originalStateIndex = originalAlignment.getState(i,j);

                // get the certain nucleotide index for each site
                int stateIndex = getAmbiguousStateIndex(originalStateIndex);

                // convert the nucleotide states into phased genotypes
                int index = getPhasedGenotypeIndex(stateIndex,stateIndex);

                // map the new alignment states
                genotypeAlignment.setState(i,j,index);
            }
        }
        return new RandomVariable<>(null, genotypeAlignment, this);
    }

    // for unit test use

    /**
     *
     * @param stateIndex state index of the nucleotide
     * @return the certain homozygous phased genotype state index
     */
    public int getAmbiguousStateIndex(int stateIndex) {
        if (stateIndex >= 4){
            // get the array for the states
            int[] ambiguousState = ambiguousState(stateIndex);
            // get the Value<Integer> for the lower and upper boundary
            Value<Integer> lower = new Value<>("id", 0);
            Value<Integer> upper = new Value<>("id",ambiguousState.length-1);

            // get the random index for the integer in the array
            UniformDiscrete uniformDiscrete = new UniformDiscrete(lower, upper);
            RandomVariable<Integer> randomNumber = uniformDiscrete.sample();

            // give the stateIndex its certain state
            stateIndex = ambiguousState[randomNumber.value()];
        }
        return stateIndex;
    }

    /**
     *
     * @param stateIndex the state index of nucleotide
     * @return the array of all possible states indices of the ambiguous nucleotide states (unkown and gap states have
     * all four possible states)
     */

    private int[] ambiguousState(int stateIndex) {
        // switch the ambiguous states into canonical states (0=A, 1=C, 2=G, 3=T)
        switch (stateIndex){
            case 4:
                // 4 = A/G
                return new int[]{0, 2};
            case 5:
                // 5 = C/T
                return new int[]{1, 3};
            case 6:
                // 6 = A/C
                return new int[]{0, 1};
            case 7:
                // 7 = A/T
                return new int[]{0, 3};
            case 8:
                // 8 = C/G
                return new int[]{1, 2};
            case 9:
                // 9 = G/T
                return new int[]{2, 3};
            case 10:
                // 10 = C/G/T
                return new int[]{1, 2, 3};
            case 11:
                // 11 = A/G/T
                return new int[]{0, 2, 3};
            case 12:
                // 12 = A/C/T
                return new int[]{0, 1, 3};
            case 13:
                // 13 = A/C/G
                return new int[]{0, 1, 2};
            case 14, 16, 15:
                // 14 = unkown base (N) = A/C/G/T
                // 15 = unkown base (?) = A/C/G/T
                // 16 = gap (-) = A/C/G/T
                return new int[]{0, 1, 2, 3};
            default:
                throw new IllegalArgumentException("Unexpected state: " + stateIndex);
        }
    }

    @Override
    public Map<String, Value> getParams() {
        return new TreeMap<>(){{
            put(ReaderConst.ALIGNMENT,alignmentValue);
        }};
    }
}
