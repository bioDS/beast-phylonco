package phylonco.beast.evolution.likelihood;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.likelihood.BeerLikelihoodCore;
import beast.base.evolution.likelihood.LikelihoodCore;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ForkJoinTask;

/**
 * {@link ReadCountTreeLikelihood} with the per-site and per-cell work of one likelihood evaluation
 * split across threads: the transition matrices of the changed branches (one task per branch), the
 * pruning at each internal node and the root sum (blocks of sites), the rebuild of tip partials
 * and of the cached read-count terms (by cell), and the read-count model's log-gamma tables (by cell
 * or read count).
 *
 * <p>Every site and every cell is computed with the same arithmetic as in ReadCountTreeLikelihood,
 * and every sum over sites or cells is taken serially in the same order, so the log likelihood does
 * not depend on the number of threads and equals ReadCountTreeLikelihood's exactly.</p>
 */
@Description("ReadCountTreeLikelihood with the work of each evaluation split across threads")
public class ThreadedReadCountTreeLikelihood extends ReadCountTreeLikelihood {

    final public Input<Integer> threadsInput = new Input<>("threads",
            "number of threads for one likelihood evaluation (default 1: the same as ReadCountTreeLikelihood)", 1);

    private int nThreads = 1;
    private ForkJoinPool pool;

    // transition matrices computed before the traversal, by node number
    private double[][] matrixBuffer;
    private double[] branchTimes;
    private boolean[] matrixNeeded;

    @Override
    public void initAndValidate() {
        nThreads = threadsInput.get();
        if (nThreads < 1) {
            throw new IllegalArgumentException("threads must be at least 1, got " + nThreads);
        }
        // the pool must exist before super.initAndValidate builds the likelihood core and tip partials
        pool = nThreads > 1 ? new ForkJoinPool(nThreads) : null;
        super.initAndValidate();
        if (pool != null && m_siteModel.getCategoryCount() != 1) {
            throw new IllegalArgumentException("ThreadedReadCountTreeLikelihood with threads > 1 is implemented for a "
                    + "single rate category only");
        }
        if (pool != null) {
            // the read-count model's table rebuilds (coverage by cell, Dirichlet by read count) on the pool
            readCountModelInput.get().setRangeRunner((n, body) -> forBlocks(n, body::accept));
        }
        int nNodes = treeInput.get().getNodeCount();
        int S = substitutionModel.getStateCount();
        matrixBuffer = new double[nNodes][S * S];
        branchTimes = new double[nNodes];
        matrixNeeded = new boolean[nNodes];
    }

    /** Runs body over [0, n) split into contiguous blocks, one per thread, and waits for all of them. */
    private void forBlocks(int n, BlockBody body) {
        int k = pool == null ? 1 : Math.min(nThreads, n);
        if (k <= 1) {
            body.run(0, n);
            return;
        }
        List<ForkJoinTask<?>> tasks = new ArrayList<>(k);
        for (int b = 0; b < k; b++) {
            final int from = (int) ((long) b * n / k), to = (int) ((long) (b + 1) * n / k);
            tasks.add(pool.submit(() -> body.run(from, to)));
        }
        for (ForkJoinTask<?> t : tasks) t.join();
    }

    @FunctionalInterface
    private interface BlockBody {
        void run(int from, int to);
    }

    // ------------------------------------------------------------------ tips

    @Override
    protected void refillTerms(List<Node> leaves) {
        forBlocks(leaves.size(), (from, to) -> super.refillTerms(leaves.subList(from, to)));
    }

    @Override
    protected double[][] computeLeafPartials(List<Node> leaves) {
        double[][] out = new double[leaves.size()][];
        forBlocks(leaves.size(), (from, to) -> {
            for (int k = from; k < to; k++) out[k] = computeLeafPartials(leaves.get(k));
        });
        return out;
    }

    // ------------------------------------------------------------------ pruning and root sum

    @Override
    protected LikelihoodCore createLikelihoodCore(int stateCount) {
        return new BlockedCore(stateCount);
    }

    /** BeerLikelihoodCore with the pruning of partials and the root sum split into blocks of sites. */
    private final class BlockedCore extends BeerLikelihoodCore {

        BlockedCore(int stateCount) {
            super(stateCount);
        }

        @Override
        public void calculatePartials(int nodeIndex1, int nodeIndex2, int nodeIndex3) {
            if (states[nodeIndex1] != null || states[nodeIndex2] != null || nrOfMatrices != 1) {
                super.calculatePartials(nodeIndex1, nodeIndex2, nodeIndex3);
                return;
            }
            final double[] p1 = partials[currentPartialsIndex[nodeIndex1]][nodeIndex1];
            final double[] m1 = matrices[currentMatrixIndex[nodeIndex1]][nodeIndex1];
            final double[] p2 = partials[currentPartialsIndex[nodeIndex2]][nodeIndex2];
            final double[] m2 = matrices[currentMatrixIndex[nodeIndex2]][nodeIndex2];
            final double[] p3 = partials[currentPartialsIndex[nodeIndex3]][nodeIndex3];
            forBlocks(nrOfPatterns, (from, to) -> prune(p1, m1, p2, m2, p3, from, to));
            if (useScaling) {
                scalePartials(nodeIndex3);
            }
        }

        /** BeerLikelihoodCore.calculatePartialsPartialsPruning (one matrix) for sites from..to-1. */
        private void prune(double[] partials1, double[] matrices1, double[] partials2, double[] matrices2,
                           double[] partials3, int from, int to) {
            final int S = nrOfStates;
            int u = from * S;
            int v = from * S;
            for (int k = from; k < to; k++) {
                int w = 0;
                for (int i = 0; i < S; i++) {
                    double sum1 = 0.0, sum2 = 0.0;
                    for (int j = 0; j < S; j++) {
                        sum1 += matrices1[w] * partials1[v + j];
                        sum2 += matrices2[w] * partials2[v + j];
                        w++;
                    }
                    partials3[u] = sum1 * sum2;
                    u++;
                }
                v += S;
            }
        }

        @Override
        public void calculateLogLikelihoods(double[] partials, double[] frequencies, double[] outLogLikelihoods) {
            forBlocks(outLogLikelihoods.length, (from, to) ->
                    rootLogLikelihoods(this, partials, frequencies, outLogLikelihoods, from, to));
        }
    }

    // ------------------------------------------------------------------ transition matrices

    /**
     * TreeLikelihood.traverse with the transition matrices of all changed branches computed first, in
     * parallel (the same getTransitionProbabilities calls), then the usual post-order traversal.
     */
    @Override
    protected int traverse(Node node) {
        if (pool == null || !node.isRoot()) {
            return super.traverse(node);
        }
        TreeInterface tree = treeInput.get();
        List<Node> changed = new ArrayList<>();
        for (Node n : tree.getNodesAsArray()) {
            if (n.isRoot()) continue;
            int nr = n.getNr();
            double branchRate = branchRateModel.getRateForBranch(n);
            branchTimes[nr] = n.getLength() * branchRate;
            int update = n.isDirty() | hasDirt;
            matrixNeeded[nr] = update != Tree.IS_CLEAN || branchTimes[nr] != m_branchLengths[nr];
            if (matrixNeeded[nr]) changed.add(n);
        }
        final double[] rates = new double[tree.getNodeCount()];
        for (Node n : changed) {
            rates[n.getNr()] = m_siteModel.getRateForCategory(0, n) * branchRateModel.getRateForBranch(n);
        }
        forBlocks(changed.size(), (from, to) -> {
            for (int k = from; k < to; k++) {
                Node n = changed.get(k);
                substitutionModel.getTransitionProbabilities(n, n.getParent().getHeight(), n.getHeight(),
                        rates[n.getNr()], matrixBuffer[n.getNr()]);
            }
        });
        return traverseWithMatrices(node);
    }

    private int traverseWithMatrices(Node node) {
        int update = node.isDirty() | hasDirt;
        final int nodeIndex = node.getNr();
        if (!node.isRoot() && matrixNeeded[nodeIndex]) {
            m_branchLengths[nodeIndex] = branchTimes[nodeIndex];
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            likelihoodCore.setNodeMatrix(nodeIndex, 0, matrixBuffer[nodeIndex]);
            update |= Tree.IS_DIRTY;
        }
        if (!node.isLeaf()) {
            final Node child1 = node.getLeft();
            final int update1 = traverseWithMatrices(child1);
            final Node child2 = node.getRight();
            final int update2 = traverseWithMatrices(child2);
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {
                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();
                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
                }
                if (m_siteModel.integrateAcrossCategories()) {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                }
                if (node.isRoot()) {
                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);
                    if (getConstantPattern() != null) {
                        proportionInvariant = m_siteModel.getProportionInvariant();
                        for (final int i : getConstantPattern()) {
                            m_fRootPartials[i] += proportionInvariant;
                        }
                    }
                    double[] rootFrequencies = substitutionModel.getFrequencies();
                    if (rootFrequenciesInput.get() != null) {
                        rootFrequencies = rootFrequenciesInput.get().getFreqs();
                    }
                    likelihoodCore.calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);
                }
            }
        }
        return update;
    }
}
