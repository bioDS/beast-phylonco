package phylonco.beast.evolution.readcountmodel;

import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
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
    public Input<Boolean> sampleAllSequencesInput = new Input<>("sampleAllSequences", "if true, sample all Sequences in one proposal; if false (default), sample one random sequence", false);
    private MATreeLikelihood maTreeLikelihood;
    private MutableAlignment mutableAlignment;
    public LikelihoodReadCountModel likelihoodReadCountModel;
    private boolean sampleAllSequences;
    private int numStates;
    private int numSites;
    private int numTaxa;
    private List<int[]> statesSequences;
    // Alignment column index -> tree node nr. Built once at init so that
    // the random Gibbs index (in [0, numTaxa)) can be interpreted as an
    // alignment column and translated to a tree node nr when needed.
    // tree-node-nr is the right index for MATreeLikelihood probes (likelihoodCore
    // is keyed by tree node); alignment-column is the right index for
    // mutableAlignment.{get,set}SiteValuesByTaxon and likelihoodReadCountModel
    // .sequenceLogLikelihood. Conflating them silently corrupts the alignment
    // on datasets where tree-leaf-order != alignment-input-order (e.g. CRC09).
    private int[] alignIdxToTreeNodeNr;




    @Override
    public void initAndValidate() {
        mutableAlignment = mutableAlignmentInput.get();
        maTreeLikelihood = maTreeLikelihoodInput.get();
        likelihoodReadCountModel = likelihoodReadCountModelInput.get();
        sampleAllSequences = sampleAllSequencesInput.get();
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

        // Cache alignment-column -> tree-node-nr mapping (via taxon name).
        TreeInterface tree = maTreeLikelihood.treeInput.get();
        alignIdxToTreeNodeNr = new int[numTaxa];
        for (int a = 0; a < numTaxa; a++) {
            String taxonName = mutableAlignment.getTaxaNames().get(a);
            int found = -1;
            for (Node leaf : tree.getExternalNodes()) {
                if (taxonName.equals(leaf.getID())) {
                    found = leaf.getNr();
                    break;
                }
            }
            if (found < 0) {
                throw new RuntimeException("GibbsSequenceOperator: alignment taxon '" +
                        taxonName + "' not found in tree leaves");
            }
            alignIdxToTreeNodeNr[a] = found;
        }
    }



    @Override
    public double proposal() {
        if(sampleAllSequences){
            int[] randomTaxaOrder = generateRandomOrder(numTaxa);
            for(int i = 0; i < numTaxa; i++){
                int kAlignIdx = randomTaxaOrder[i];
                int kNodeNr = alignIdxToTreeNodeNr[kAlignIdx];
                int[] newSeq = sampleSequence(kNodeNr, kAlignIdx);
                mutableAlignment.setSiteValuesByTaxon(kAlignIdx, newSeq);
                maTreeLikelihood.getLogProbsForStateSequence(kNodeNr, newSeq);
            }
        } else {
            int kAlignIdx = Randomizer.nextInt(numTaxa);
            int kNodeNr = alignIdxToTreeNodeNr[kAlignIdx];
            int[] newSeq = sampleSequence(kNodeNr, kAlignIdx);
            mutableAlignment.setSiteValuesByTaxon(kAlignIdx, newSeq);
        }
        return Double.POSITIVE_INFINITY;
    }

    private int[] sampleSequence(int kNodeNr, int kAlignIdx) {
        double[][] stateLogProbabilities = new double[numStates][numSites];
        double[][] readCountLogLikelihoods = new double[numStates][numSites];
        int[] newSeq = new int[numSites];
        double[] stateProbabilities;

        for(int i = 0; i < numStates; i++){
            // tree partials indexed by tree node nr
            stateLogProbabilities[i] = maTreeLikelihood.getLogProbsForStateSequence(kNodeNr, statesSequences.get(i));
            // read counts indexed by alignment column
            readCountLogLikelihoods[i] = likelihoodReadCountModel.sequenceLogLikelihood(kAlignIdx, statesSequences.get(i));
        }

        for (int i = 0; i < numSites; i++) {
            // multiple read count likelihoods and statelogLikelihoods
            // * for this site * and normalise and sample from and set sequence at site i
            double[] logProbs = new double[numStates];
            for(int j = 0; j < numStates; j++){
                logProbs[j] = stateLogProbabilities[j][i] + readCountLogLikelihoods[j][i];
            }
            stateProbabilities = normalizeLogProbs(logProbs);
            newSeq[i] = sampleStateFromProbabilities(stateProbabilities);
        }
        return newSeq;
    }

    private int sampleStateFromProbabilities(double[] probabilities) {
        double rand = Randomizer.nextDouble();
        double cumulative = 0.0;
        for (int i = 0; i < probabilities.length; i++) {
            cumulative += probabilities[i];
            if (rand < cumulative) return i;
        }
        return sampleStateFromProbabilities(probabilities);
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

    public static int[] generateRandomOrder(int n) {
        int[] result = new int[n];
        for (int i = 0; i < n; i++) {
            result[i] = i;
        }

        for (int i = n - 1; i > 0; i--) {
            int j = Randomizer.nextInt(i+1);
            int temp = result[i];
            result[i] = result[j];
            result[j] = temp;
        }

        return result;
    }
}
