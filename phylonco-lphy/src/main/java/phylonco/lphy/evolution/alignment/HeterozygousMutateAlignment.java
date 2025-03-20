package phylonco.lphy.evolution.alignment;

import jebl.evolution.sequences.Nucleotides;
import lphy.base.distribution.ParametricDistribution;
import lphy.base.evolution.SNPSampler;
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

import java.util.*;

import static lphy.base.evolution.SNPSampler.getAmbiguousStateIndex;
import static lphy.base.evolution.SNPSampler.sampleRandomNumber;
import static phylonco.lphy.evolution.datatype.PhasedGenotype.getNucleotideIndex;
import static phylonco.lphy.evolution.datatype.PhasedGenotype.getPhasedGenotypeIndex;

public class HeterozygousMutateAlignment extends ParametricDistribution<Alignment> {
    private Value<Alignment> alignment;
    private Value<Integer> n;
    private Value<Integer[]> positions;
    private static final String numberName = "n";
    private static final String positionName = "positions";

    int altIndex;
    int refIndex;

    public HeterozygousMutateAlignment(
            @ParameterInfo(name = ReaderConst.ALIGNMENT, description = "the input alignment") Value<Alignment> alignment,
            @ParameterInfo(name = numberName, description = "the number of heterozygous sites") Value<Integer> n,
            @ParameterInfo(name = positionName, description = "positions that heterozygous sites occur, should be fewer than the number", optional = true) Value<Integer[]> positions){
        if (alignment == null) {
            throw new IllegalArgumentException("Alignment cannot be null");
        }
        if (n == null) {
            throw new IllegalArgumentException("Number of sites cannot be null");
        }
        if (positions != null && positions.value().length > n.value()){
            throw new IllegalArgumentException("Number of sites (n) should be equal or greater than the number of positions (positions)!");
        }
        //TODO: handle  length < 1000
        this.alignment = alignment;
        this.n = n;
        this.positions = positions;
    }

    @Override
    protected void constructDistribution(RandomGenerator random) {

    }

    @GeneratorInfo(name = "Heterozygous", description = "Sample n sites in the given alignment and convert the sites to heterozygous ref/alt. If the site is homozygous, take any allele as ref then randomly sample an alt." +
            "If the site is heterozygous, take any allele as ref and randomly sample an alt. If alignment is haploid, then convert the sites to heterozygoues and others to homozygous. ")
    @Override
    public RandomVariable<Alignment> sample() {

        // get all params we need
        Alignment alignment = getAlignment().value();
        int length = alignment.nchar();
        int n = getNumber().value();
        Integer[] positions = new Integer[]{};
        if ( getPositions() != null) {
           positions = getPositions().value();
        }

        // obtain the taxa names
        List<String> taxa = new ArrayList<>();
        for (int s = 0; s < alignment.ntaxa(); s++) {
            taxa.add(alignment.getTaxonName(s));
        }
        String[] taxaNames = taxa.toArray(new String[0]);

        // initialise the output alignment
        Alignment outAlignment = new SimpleAlignment(Taxa.createTaxa(taxaNames),
                alignment.nchar(), PhasedGenotype.INSTANCE);

        // set the alignment
        for (int i = 0; i < outAlignment.ntaxa(); i++) {
            // construct positions to convert
            List<Integer> positionSet  = List.of(positions);
            if (n != positions.length) {
                int extraNumber = n - positions.length;
                positionSet = samplePositions(positions, extraNumber, length);
            }

            for (int j = 0; j < outAlignment.nchar(); j++) {
                // get the nucleotide index of each site
                int inputIndex = alignment.getState(i, j);
                int outputIndex;
                // check whether the site is in positionSet
                // if not in the positionSet
                if (! positionSet. contains(j)){
                    // convert to homozygous if its nucleotide
                    if (alignment.getSequenceTypeStr().equals(Nucleotides.NAME)){
                        inputIndex = getAmbiguousStateIndex(inputIndex);
                        outputIndex = getPhasedGenotypeIndex(inputIndex, inputIndex);
                    } else {
                        inputIndex = getGT16AmbiguousStateIndex(inputIndex);
                        outputIndex = inputIndex;
                    }
                    outAlignment.setState(i, j, outputIndex);
                    // if in the positionSet
                } else{
                    // if haploid
                    // consider the state as ref and sample an alt
                    if (alignment.getSequenceTypeStr().equals(Nucleotides.NAME)){
                        inputIndex = getAmbiguousStateIndex(inputIndex);
                        refIndex = inputIndex;
                        altIndex = getRandomCanonicalState(new int[]{inputIndex});
                    } else {
                        inputIndex = getGT16AmbiguousStateIndex(inputIndex);
                        int[] parentIndices = getNucleotideIndex(inputIndex);
                        // if homozygous
                        if (parentIndices[0] == parentIndices[1]) {
                            refIndex = parentIndices[0];
                            altIndex = getRandomCanonicalState(new int[]{parentIndices[0]});
                        } else {
                            // if heterozygous
                            refIndex = parentIndices[sampleRandomNumber(0, parentIndices.length - 1)];
                            altIndex = getRandomCanonicalState(parentIndices);
                        }
                    }
                    outputIndex = getPhasedGenotypeIndex(refIndex, altIndex);
                    outAlignment.setState(i, j, outputIndex);
                }
            }
        }
        return new RandomVariable<>(null, outAlignment, this);
    }

    private int getGT16AmbiguousStateIndex(int inputIndex) {
        if (inputIndex>=0 && inputIndex<=15){
            return inputIndex;
        } else {
            return sampleRandomNumber(0,15);
        }
    }

    public static int getRandomCanonicalState(int[] refIndex) {
        List<Integer> stateIndices = new ArrayList<>(Arrays.asList(0, 1, 2, 3));
        for (int index : refIndex) {
            stateIndices.remove(Integer.valueOf(index));
        }

        int randomIndex = sampleRandomNumber(0, stateIndices.size() - 1);
        return stateIndices.get(randomIndex);
    }

    public static List<Integer> samplePositions(Integer[] positions, int extraNumber, int alignmentLength) {
        List<Integer> positionSet = new ArrayList<>(List.of(positions));
        for (int i = 0; i < extraNumber; i++){
            int newPosition = sampleRandomNumber(0, alignmentLength - 1);
            while (positionSet.contains(newPosition)){
                newPosition = sampleRandomNumber(0, alignmentLength - 1);
            }
            positionSet.add(newPosition);
        }
        Collections.sort(positionSet);
        return positionSet;
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> params = new TreeMap<>();
        if (alignment != null) params.put(ReaderConst.ALIGNMENT, alignment);
        if (n != null) params.put(numberName, n);
        if (positions != null) params.put(positionName, positions);
       return params;
    }

    public void setParam(String paramName, Value value){
        if (paramName.equals(ReaderConst.ALIGNMENT)) alignment = value;
        else if (paramName.equals(numberName)) n = value;
        else if (paramName.equals(positionName)) positions = value;
    }

    public Value<Alignment> getAlignment(){
        return getParams().get(ReaderConst.ALIGNMENT);
    }

    public Value<Integer> getNumber(){
        return getParams().get(numberName);
    }

    public Value<Integer[]> getPositions(){
        return getParams().get(positionName);
    }
}
