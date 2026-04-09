package phylonco.beast.evolution.readcountmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import mutablealignment.MATreeLikelihood;
import mutablealignment.MutableAlignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Compound operator: narrow exchange followed by Gibbs resampling of one
 * affected leaf taxon's genotype sequence.
 *
 * <p>When tree topology and genotypes are strongly coupled, proposing a new
 * tree alone is unlikely to be accepted because the current genotypes are
 * tuned to the old tree. This operator jointly proposes a topology change
 * and resamples the genotypes of one affected taxon, giving the new tree
 * a fair evaluation.</p>
 *
 * <h3>Acceptance ratio</h3>
 *
 * <p>For a narrow exchange swapping nodes i and uncle, with affected leaf
 * taxon k whose sequence is Gibbs-resampled, the Hastings ratio returned
 * by this operator is:</p>
 *
 * <pre>
 *   logHR = logHR_exchange + logQ_reverse - logQ_forward
 * </pre>
 *
 * <p>where logHR_exchange = log(validGP / validGPafter) is the standard
 * narrow-exchange Hastings ratio, and</p>
 *
 * <pre>
 *   logQ_forward = sum_s log q(a'_k[s] | T', A_{-k}, s)   (new seq, new tree)
 *   logQ_reverse = sum_s log q(a_k[s]  | T,  A_{-k}, s)   (old seq, old tree)
 * </pre>
 *
 * <p>with q(g | T, A_{-k}, s) = P(R_k[s]|g) * L_s(A_{-k}, g | T) / Z_s
 * being the normalised Gibbs probability at site s.</p>
 *
 * <p>The target-distribution terms cancel with the Gibbs proposal terms,
 * so the full acceptance ratio simplifies to:</p>
 *
 * <pre>
 *   alpha = P(T')/P(T) * (validGP/validGPafter) * prod_s Z_s(T',A_{-k}) / Z_s(T,A_{-k})
 * </pre>
 */
@Description("Narrow exchange + Gibbs resample of one affected leaf taxon. " +
        "Improves mixing when tree topology and genotypes are strongly coupled.")
public class ExchangeGibbsOperator extends TreeOperator {

    public Input<MutableAlignment> mutableAlignmentInput = new Input<>(
            "mutableAlignment",
            "mutable alignment containing genotype states",
            Input.Validate.REQUIRED);

    public Input<MATreeLikelihood> maTreeLikelihoodInput = new Input<>(
            "maTreeLikelihood",
            "tree likelihood for mutable alignment",
            Input.Validate.REQUIRED);

    public Input<LikelihoodReadCountModel> likelihoodReadCountModelInput = new Input<>(
            "likelihoodReadCountModel",
            "read count likelihood model",
            Input.Validate.REQUIRED);

    private MutableAlignment mutableAlignment;
    private MATreeLikelihood maTreeLikelihood;
    private LikelihoodReadCountModel likelihoodReadCountModel;

    private int numStates;
    private int numSites;
    private List<int[]> allGSequences;

    @Override
    public void initAndValidate() {
        mutableAlignment = mutableAlignmentInput.get();
        maTreeLikelihood = maTreeLikelihoodInput.get();
        likelihoodReadCountModel = likelihoodReadCountModelInput.get();
        numStates = mutableAlignment.getDataType().getStateCount();
        numSites = mutableAlignment.getSiteCount();

        // Precompute constant sequences: allGSequences.get(g) = [g, g, ..., g]
        allGSequences = new ArrayList<>(numStates);
        for (int g = 0; g < numStates; g++) {
            int[] seq = new int[numSites];
            Arrays.fill(seq, g);
            allGSequences.add(seq);
        }
    }

    @Override
    public double proposal() {
        final Tree tree = (Tree) InputUtil.get(treeInput, this);
        final int internalNodes = tree.getInternalNodeCount();

        if (internalNodes <= 1) {
            return Double.NEGATIVE_INFINITY;
        }

        // ---- Pick grandparent (internal node with at least one non-leaf child) ----
        Node grandParent = tree.getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
        while (grandParent.getLeft().isLeaf() && grandParent.getRight().isLeaf()) {
            grandParent = tree.getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
        }

        // parentIndex = higher child, uncle = lower child
        Node parentIndex = grandParent.getLeft();
        Node uncle = grandParent.getRight();
        if (parentIndex.getHeight() < uncle.getHeight()) {
            parentIndex = grandParent.getRight();
            uncle = grandParent.getLeft();
        }

        if (parentIndex.isLeaf()) {
            // tree with dated tips
            return Double.NEGATIVE_INFINITY;
        }

        // ---- Count valid grandparents before exchange ----
        int validGP = 0;
        for (int n = internalNodes + 1; n < 1 + 2 * internalNodes; n++) {
            validGP += isg(tree.getNode(n));
        }
        final int c2 = sisg(parentIndex) + sisg(uncle);

        // Pick random child of parentIndex to swap with uncle
        final Node i = Randomizer.nextBoolean() ? parentIndex.getLeft() : parentIndex.getRight();

        // ---- Collect affected leaves from both swapped subtrees ----
        List<Node> iLeaves = collectLeaves(i);
        List<Node> uncleLeaves = collectLeaves(uncle);
        List<Node> allLeaves = new ArrayList<>(iLeaves);
        allLeaves.addAll(uncleLeaves);

        if (allLeaves.isEmpty()) {
            return Double.NEGATIVE_INFINITY;
        }

        // Pick one affected leaf taxon k
        Node leafNode = allLeaves.get(Randomizer.nextInt(allLeaves.size()));
        int k = leafNode.getNr();
        boolean kFromI = iLeaves.contains(leafNode);

        // ---- Compute logQ_reverse (under current tree T, before exchange) ----
        int[] currentSeq = mutableAlignment.getSiteValuesByTaxon(k);
        double logQ_reverse = computeLogGibbsProb(k, currentSeq);

        // ---- Perform narrow exchange: T -> T' ----
        exchangeNodes(i, uncle, parentIndex, grandParent);
        final int validGPafter = validGP - c2 + sisg(parentIndex) + sisg(uncle);
        double logHR_exchange = Math.log((float) validGP / validGPafter);

        // ---- Fix cached partials after topology change ----
        // After the exchange, parentIndex's children changed (uncle moved in,
        // i moved out). If k came from i's subtree, the forward Gibbs
        // computation walks from k to root WITHOUT passing through parentIndex,
        // so parentIndex's stale partials would corrupt the result.
        // Fix: walk a leaf under uncle (now a child of parentIndex) to the root,
        // which recalculates parentIndex's partials from its new children.
        if (kFromI) {
            Node fixLeaf = uncleLeaves.get(0);
            int fixNr = fixLeaf.getNr();
            maTreeLikelihood.getLogProbsForStateSequence(
                    fixNr, mutableAlignment.getSiteValuesByTaxon(fixNr));
        }

        // ---- Gibbs resample + compute logQ_forward (under new tree T') ----
        double[][] treeLogProb = new double[numStates][];
        double[][] rcLogLik = new double[numStates][];
        for (int g = 0; g < numStates; g++) {
            treeLogProb[g] = maTreeLikelihood.getLogProbsForStateSequence(k, allGSequences.get(g));
            rcLogLik[g] = likelihoodReadCountModel.sequenceLogLikelihood(k, allGSequences.get(g));
        }

        int[] newSeq = new int[numSites];
        double logQ_forward = 0.0;
        for (int s = 0; s < numSites; s++) {
            double[] logP = new double[numStates];
            for (int g = 0; g < numStates; g++) {
                logP[g] = treeLogProb[g][s] + rcLogLik[g][s];
            }
            double[] probs = normalizeLogProbs(logP);
            newSeq[s] = sampleFromProbabilities(probs);
            logQ_forward += Math.log(probs[newSeq[s]]);
        }
        mutableAlignment.setSiteValuesByTaxon(k, newSeq);

        // ---- Return compound Hastings ratio ----
        return logHR_exchange + logQ_reverse - logQ_forward;
    }

    /**
     * Compute log q_G(seq | T, A_{-k}) = sum_s log q(seq[s] | T, A_{-k}, s).
     *
     * Uses the current tree state in maTreeLikelihood (call before or after
     * exchange as appropriate).
     */
    private double computeLogGibbsProb(int k, int[] seq) {
        double[][] treeLogProb = new double[numStates][];
        double[][] rcLogLik = new double[numStates][];
        for (int g = 0; g < numStates; g++) {
            treeLogProb[g] = maTreeLikelihood.getLogProbsForStateSequence(k, allGSequences.get(g));
            rcLogLik[g] = likelihoodReadCountModel.sequenceLogLikelihood(k, allGSequences.get(g));
        }

        double logProb = 0.0;
        for (int s = 0; s < numSites; s++) {
            double[] logP = new double[numStates];
            for (int g = 0; g < numStates; g++) {
                logP[g] = treeLogProb[g][s] + rcLogLik[g][s];
            }
            double logZ = logSumExp(logP);
            logProb += logP[seq[s]] - logZ;
        }
        return logProb;
    }

    /**
     * Swap subtrees rooted at nodes i and j, whose parents are p and jP.
     * Same logic as Exchange.exchangeNodes.
     */
    private void exchangeNodes(Node i, Node j, Node p, Node jP) {
        replace(p, i, j);
        replace(jP, j, i);
    }

    /**
     * Collect all leaf nodes under the given node (including the node itself
     * if it is a leaf).
     */
    private List<Node> collectLeaves(Node node) {
        List<Node> leaves = new ArrayList<>();
        if (node.isLeaf()) {
            leaves.add(node);
        } else {
            // Node.getAllLeafNodes(List) adds leaf descendants recursively
            node.getAllLeafNodes(leaves);
        }
        return leaves;
    }

    /** 1 if node has at least one non-leaf child, 0 otherwise. */
    private int isg(Node n) {
        return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
    }

    /** isg for a node that might itself be a leaf (returns 0 for leaves). */
    private int sisg(Node n) {
        return n.isLeaf() ? 0 : isg(n);
    }

    private double logSumExp(double[] logP) {
        double max = Double.NEGATIVE_INFINITY;
        for (double v : logP) {
            if (v > max) max = v;
        }
        double sum = 0.0;
        for (double v : logP) {
            sum += Math.exp(v - max);
        }
        return max + Math.log(sum);
    }

    private double[] normalizeLogProbs(double[] logP) {
        double max = Double.NEGATIVE_INFINITY;
        for (double v : logP) {
            if (v > max) max = v;
        }
        double[] probs = new double[logP.length];
        double sum = 0.0;
        for (int j = 0; j < logP.length; j++) {
            probs[j] = Math.exp(logP[j] - max);
            sum += probs[j];
        }
        for (int j = 0; j < probs.length; j++) {
            probs[j] /= sum;
        }
        return probs;
    }

    private int sampleFromProbabilities(double[] probs) {
        double rand = Randomizer.nextDouble();
        double cumulative = 0.0;
        for (int j = 0; j < probs.length; j++) {
            cumulative += probs[j];
            if (rand < cumulative) return j;
        }
        return probs.length - 1;
    }
}
