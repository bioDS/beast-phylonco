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



public class GibbsAlignmentOperator extends Operator {

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
    // Alignment column index -> tree node nr (via taxon name).
    // Needed because likelihoodCore probes are keyed by tree node, while
    // alignment access and read-count likelihoods are keyed by alignment column.
    private int[] alignIdxToTreeNodeNr;




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
                throw new RuntimeException("GibbsAlignmentOperator: alignment taxon '" +
                        taxonName + "' not found in tree leaves");
            }
            alignIdxToTreeNodeNr[a] = found;
        }
    }



    @Override
    public double proposal() {
        double[][] stateLogProbabilities = new double[numStates][numSites];
        double[][] readCountLogLikelihoods = new double[numStates][numSites];
        int[] newSeq = new int[numSites];
        double[] stateProbabilities;
        int[] randomTaxaOrder = generateRandomOrder(numTaxa);
        for (int k = 0; k < numTaxa; k++) {
            int kAlignIdx = randomTaxaOrder[k];
            int kNodeNr = alignIdxToTreeNodeNr[kAlignIdx];
            for (int i = 0; i < numStates; i++) {
                // tree partials indexed by tree node nr
                stateLogProbabilities[i] = maTreeLikelihood.getLogProbsForStateSequence(kNodeNr, statesSequences.get(i));
                // read counts indexed by alignment column
                readCountLogLikelihoods[i] = likelihoodReadCountModel.sequenceLogLikelihood(kAlignIdx, statesSequences.get(i));
            }

            for (int i = 0; i < numSites; i++) {
                // multiple read count likelihoods and statelogLikelihoods
                // * for this site * and normalise and sample from and set sequence at site i
                double[] logProbs = new double[numStates];
                for (int j = 0; j < numStates; j++) {
                    logProbs[j] = stateLogProbabilities[j][i] + readCountLogLikelihoods[j][i];
                }
                stateProbabilities = normalizeLogProbs(logProbs);
                newSeq[i] = sampleFromProbabilities(stateProbabilities);
            }


            mutableAlignment.setSiteValuesByTaxon(kAlignIdx, newSeq);
            // needed to update partial log likelihoods based on newly sampled sequence before
            // moving to next taxa. Don't need to use the result
            maTreeLikelihood.getLogProbsForStateSequence(kNodeNr, newSeq);
        }
        return Double.POSITIVE_INFINITY;
    }

    private int sampleFromProbabilities(double[] probabilities) {
        double rand = Randomizer.nextDouble();
        double cumulative = 0.0;
        for (int i = 0; i < probabilities.length; i++) {
            cumulative += probabilities[i];
            if (rand < cumulative) return i;
        }
        throw new RuntimeException("Should never get here");
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
