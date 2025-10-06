package phylonco.beast.evolution.readcountmodel;

import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.util.Randomizer;
import mutablealignment.MATreeLikelihood;
import mutablealignment.MutableAlignment;

import java.util.ArrayList;
import java.util.List;



public class GibbsSequenceOperator extends Operator {

    public Input<MutableAlignment> mutableAlignmentInput = new Input<>("mutableAlignment", "mutable alignment");
    public Input<MATreeLikelihood> maTreeLikelihoodInput  = new Input<>("maTreeLikelihood", "likelihood of mutable alignment tree");
    public Input<LikelihoodReadCountModel> likelihoodReadCountModelInput = new Input<>("likelihoodReadCountModel", "ikelihood of Read Count Model");
    private MATreeLikelihood maTreeLikelihood;
    private MutableAlignment mutableAlignment;
    public LikelihoodReadCountModel likelihoodReadCountModel;
    private int numStates;
    private int numSites;
    private int numTaxa;
    private List<int[]> statesSequences;




    @Override
    public void initAndValidate() {
        mutableAlignment = mutableAlignmentInput.get();
        maTreeLikelihood = maTreeLikelihoodInput.get();
        likelihoodReadCountModel = likelihoodReadCountModelInput.get();
        numStates = mutableAlignment.getDataType().getStateCount();
        numSites = mutableAlignment.getSiteCount();
        numTaxa = mutableAlignment.getTaxonCount();
        statesSequences = new ArrayList<>(numStates);
        for(int i = 0; i < numStates; i++){
            int[] stateSequence = new int[numSites];
            for(int j = 0; j < numSites; j++){
                stateSequence[j] = i;
            }
            statesSequences.add(stateSequence);
        }

    }



    @Override
    public double proposal() {
        int taxon = Randomizer.nextInt(numTaxa);
        double[][] stateLogProbabilities = new double[numStates][numSites];
        double[][] readCountLogLikelihoods = new double[numStates][numSites];
        int[] newSeq = new int[numSites];
        double[] stateProbabilities;

        for(int i = 0; i < numStates; i++){
            stateLogProbabilities[i] = maTreeLikelihood.getLogProbsForStateSequence(taxon, statesSequences.get(i));

            // get new read count likelihood for taxon
            readCountLogLikelihoods[i] = likelihoodReadCountModel.sequenceLogLikelihood(taxon,statesSequences.get(i));
        }

        for (int i = 0; i < numSites; i++) {
            // multiple read count likelihoods and statelogLikelihoods
            // * for this site * and normalise and sample from and set sequence at site i
            double[] logProbs = new double[numStates];
            for(int j = 0; j < numStates; j++){
                logProbs[j] = stateLogProbabilities[j][i] + readCountLogLikelihoods[j][i];
            }
            stateProbabilities = normalizeLogProbs(logProbs);
            newSeq[i] = sampleFromProbabilities(stateProbabilities);
        }
        mutableAlignment.setSiteValuesByTaxon(taxon, newSeq);
        return Double.POSITIVE_INFINITY;
    }

    private int sampleFromProbabilities(double[] probabilities) {
        double rand = Randomizer.nextDouble();
        double cumulative = 0.0;
        for (int i = 0; i < probabilities.length; i++) {
            cumulative += probabilities[i];
            if (rand < cumulative) return i;
        }
        return sampleFromProbabilities(probabilities);
    }

    private double[] normalizeLogProbs(double[] logProbs) {
        // Find the maximum log probability for numerical stability
        double maxLogProb = Double.NEGATIVE_INFINITY;
        for (double logP : logProbs) {
            maxLogProb = Math.max(maxLogProb, logP);
        }
        // Compute exp(logProb - max) for each element and sum them
        double[] expProbs = new double[logProbs.length];
        double sumExp = 0.0;
        for (int i = 0; i < logProbs.length; i++) {
            expProbs[i] = Math.exp(logProbs[i] - maxLogProb);
            sumExp += expProbs[i];
        }
        // Normalize to get actual probabilities
        double[] normalizedProbs = new double[logProbs.length];
        for (int i = 0; i < logProbs.length; i++) {
            normalizedProbs[i] = expProbs[i] / sumExp;
        }
        return normalizedProbs;
    }
}
