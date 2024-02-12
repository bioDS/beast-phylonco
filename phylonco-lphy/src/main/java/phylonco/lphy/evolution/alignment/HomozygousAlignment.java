package phylonco.lphy.evolution.alignment;

import jebl.evolution.sequences.Nucleotides;
import lphy.base.evolution.Taxa;
import lphy.base.evolution.alignment.AbstractAlignment;
import lphy.base.evolution.alignment.Alignment;
import lphy.base.evolution.alignment.SimpleAlignment;
import lphy.base.function.io.ReaderConst;
import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import phylonco.lphy.evolution.datatype.PhasedGenotype;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;


public class HomozygousAlignment extends DeterministicFunction<Alignment> {
    public HomozygousAlignment(@ParameterInfo(name = ReaderConst.ALIGNMENT,
            description = "the genotype alignment (homozygous) converted from diploid alignment" )
                            Value<AbstractAlignment> alignmentValue){
        if (alignmentValue == null) throw new IllegalArgumentException("The alignment can't be null!");
        setParam(ReaderConst.ALIGNMENT, alignmentValue);
    }

    // write the function
    @GeneratorInfo(name = "Homozygous", description = "Convert the haploid sequence to genotype sequence. Give the accurate " +
            "homozygous alignment when there are only canonical states in the sequence. If there are ambiguous states in the sequence," +
            "state will be chosen randomly among the possible states and convert to the homozygous alignment.")
    @Override
    public Value<Alignment> apply() {
        // get the original seq
        Alignment originalAlignment = ((Value<Alignment>) getParams().get(ReaderConst.ALIGNMENT)).value();

        // obtain the taxa names
        List<String> taxa = new ArrayList<>();
        for (int s =0; s < originalAlignment.ntaxa();s++){
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
                int stateIndex = originalAlignment.getState(i,j);

                // get the certain nucleotide index for each site
                if (stateIndex >= 4){
                    stateIndex = ambiguousState(stateIndex);
                }
                // convert the nucleotide states into phased genotypes
                int index = 4*stateIndex + stateIndex;

                // map the new alignment states
                genotypeAlignment.setState(i,j,index);
            }
        }

        return new Value <>(null, genotypeAlignment, this);
    }

    // deal with the ambiguities
    private int ambiguousState(int stateIndex) {
        // create a random object
        Random r = new Random();
        // switch the ambiguous states into canonical states (0=A, 1=C, 2=G, 3=T)
        switch (stateIndex){
            case 4:
                // 4 = A/G
                int[] R = {0,2};
                // get the random index for choosing the state
                int randomIndex_R = r.nextInt(R.length);
                // get the selected state canonical state
                int randomState_R = R[randomIndex_R];
                return randomState_R;
            case 5:
                // 5 =C/T
                int[] Y = {1,3};
                // short for choosing the random index and get the canonical state
                return Y[r.nextInt(Y.length)];
            case 6:
                // 6 = A/C
                int[] M = {0,1};
                return M[r.nextInt(M.length)];
            case 7:
                // 7 = A/T
                int[] W = {0,3};
                return W[r.nextInt(W.length)];
            case 8:
                // 8 = C/G
                int[] S = {1,2};
                return S[r.nextInt(S.length)];
            case 9:
                // 9 = G/T
                int[] K = {2,3};
                return K[r.nextInt(K.length)];
            case 10:
                // 10 = C/G/T
                int[] B = {1,2,3};
                return B[r.nextInt(B.length)];
            case 11:
                // 11 = A/G/T
                int[] D = {0,2,3};
                return D[r.nextInt(D.length)];
            case 12:
                // 12 = A/C/T
                int[] H = {0,1,3};
                return H[r.nextInt(H.length)];
            case 13:
                // 13 = A/C/G
                int[] V = {0,1,2};
                return V[r.nextInt(V.length)];
            case 14:
                // 14 = unkown base = A/C/G/T
                int[] N = {0,1,2,3};
                return N[r.nextInt(N.length)];
            case 15:
                // 15 = unkown base = A/C/G/T
                int[] unkown = {0,1,2,3};
                return unkown[r.nextInt(unkown.length)];
            case 16:
                // 16 = gap = A/C/G/T
                int[] gap = {0,1,2,3};
                return gap[r.nextInt(gap.length)];
        } throw new RuntimeException("Unexpected state: " + stateIndex);
    }
}
