package phylonco.beast.evolution.likelihood;

import beagle.Beagle;
import beagle.BeagleFactory;
import beagle.BeagleFlag;
import beagle.InstanceDetails;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;

import java.util.ArrayList;
import java.util.List;

/**
 * {@link ReadCountTreeLikelihood} with the tree part of the likelihood computed by BEAGLE (CPU, SSE
 * or GPU): the pruning of partials and their rescaling. Everything else is inherited unchanged: the
 * read-count tip partials (cached Dirichlet terms, rebuilt only for the cells that changed, split
 * across threads with {@code threads}) and genotype sampling.
 *
 * <p>The transition matrices are computed in Java exactly as ReadCountTreeLikelihood computes them
 * and handed to BEAGLE (so the genotype sampler uses the same matrices); the root
 * sum over genotypes is taken in Java from BEAGLE's root
 * partials and accumulated log scale factors. Partials, matrices and scale factors are double-buffered
 * for store/restore, as in BEAST's BeagleTreeLikelihood. Every internal node is rescaled, so values
 * agree with ReadCountTreeLikelihood to rounding (not bit for bit).</p>
 *
 * <p>The BEAGLE resource and flags follow BEAST's options (-beagle_GPU, -beagle_CPU, -beagle_SSE,
 * -beagle_double, -beagle_order, i.e. the system properties beagle.resource.order,
 * beagle.preferred.flags and beagle.required.flags); double precision is preferred unless single is
 * asked for. Needs the BEAGLE native library (java.library.path).</p>
 */
@Description("ReadCountTreeLikelihood with the pruning computed by BEAGLE (CPU or GPU)")
public class BeagleReadCountTreeLikelihood extends ThreadedReadCountTreeLikelihood {

    final public Input<Integer> beagleThreadsInput = new Input<>("beagleThreads",
            "number of CPU threads for a BEAGLE CPU instance (default 1)", 1);

    private Beagle beagle;
    private int nNodes, nTips, nInternal, S, nPatterns, nCat;
    private BufferIndexHelper partialHelper, matrixHelper, scaleHelper;
    private int cumulativeScaleIndex;
    private double[][] javaMatrices;          // by matrix buffer index: what BEAGLE holds, for the sampler
    private double[] rootPartials, logScales, matrixScratch, uniform;
    private final double[] sumLogLikelihood = new double[1];
    private int[] operations;
    private int[] scaleIndices;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        TreeInterface tree = treeInput.get();
        nNodes = tree.getNodeCount();
        nTips = tree.getLeafNodeCount();
        nInternal = nNodes - nTips;
        S = substitutionModel.getStateCount();
        nPatterns = dataInput.get().getPatternCount();
        nCat = m_siteModel.getCategoryCount();

        partialHelper = new BufferIndexHelper(nNodes, 0);
        matrixHelper = new BufferIndexHelper(nNodes, 0);
        scaleHelper = new BufferIndexHelper(nInternal, 0);
        cumulativeScaleIndex = scaleHelper.getBufferCount();
        javaMatrices = new double[matrixHelper.getBufferCount()][];
        rootPartials = new double[nCat * nPatterns * S];
        logScales = new double[nPatterns];
        matrixScratch = new double[S * S];
        operations = new int[nInternal * Beagle.OPERATION_TUPLE_SIZE];
        scaleIndices = new int[nInternal];

        beagle = createInstance();
        double[] weights = new double[nPatterns];
        java.util.Arrays.fill(weights, 1.0);
        beagle.setPatternWeights(weights);
        uniform = new double[S];
        java.util.Arrays.fill(uniform, 1.0 / S);
        beagle.setStateFrequencies(0, uniform);

        // tips: the read-count partials the base class built for the Java core; upload them once
        for (Node leaf : tree.getExternalNodes()) {
            beagle.setPartials(partialHelper.getOffsetIndex(leaf.getNr()), computeLeafPartials(leaf));
        }
        hasDirt = Tree.IS_FILTHY;
    }

    private Beagle createInstance() {
        long preference = 0, requirement = 0;
        int[] resources = null;
        List<Long> order = parseLongs(System.getProperty("beagle.resource.order"));
        if (!order.isEmpty()) {
            resources = new int[]{order.get(0).intValue(), 0};
            if (order.get(0) > 0) preference |= BeagleFlag.PROCESSOR_GPU.getMask();
        }
        List<Long> preferred = parseLongs(System.getProperty("beagle.preferred.flags"));
        if (!preferred.isEmpty()) preference = preferred.get(0);
        List<Long> required = parseLongs(System.getProperty("beagle.required.flags"));
        if (!required.isEmpty()) requirement = required.get(0);
        if ((preference & BeagleFlag.PRECISION_SINGLE.getMask()) == 0) preference |= BeagleFlag.PRECISION_DOUBLE.getMask();
        preference |= BeagleFlag.SCALING_MANUAL.getMask() | BeagleFlag.SCALERS_LOG.getMask();
        int threads = beagleThreadsInput.get();
        if (threads > 1) preference |= BeagleFlag.THREADING_CPP.getMask();

        Beagle b;
        Throwable failure = null;
        try {
            b = BeagleFactory.loadBeagleInstance(nTips, partialHelper.getBufferCount(), 0, S, nPatterns,
                    1, matrixHelper.getBufferCount(), nCat, scaleHelper.getBufferCount() + 1,
                    resources, preference, requirement);
        } catch (Throwable e) {
            b = null;
            failure = e;
        }
        if (b == null) {
            throw new IllegalArgumentException("BeagleReadCountTreeLikelihood: no BEAGLE instance (is the BEAGLE "
                    + "library on java.library.path?); use ReadCountTreeLikelihood without BEAGLE"
                    + (failure != null ? ". Cause: " + failure : "") + " [resources "
                    + java.util.Arrays.toString(resources) + ", preference flags " + preference + ", required " + requirement
                    + "]", failure);
        }
        if (threads > 1) b.setCPUThreadCount(threads);
        InstanceDetails d = b.getDetails();
        if (d != null) {
            Log.info.println("  BeagleReadCountTreeLikelihood: " + d.toString().trim());
        }
        return b;
    }

    private static List<Long> parseLongs(String s) {
        List<Long> out = new ArrayList<>();
        if (s == null || s.isBlank()) return out;
        for (String part : s.split(",")) out.add(Long.parseLong(part.trim()));
        return out;
    }

    // ------------------------------------------------------------------ engine hooks

    @Override
    protected void setLeafPartialsInCore(int nodeNr, double[] partials) {
        super.setLeafPartialsInCore(nodeNr, partials);   // the Java core too (cheap; keeps it consistent)
        if (beagle == null) return;                       // during super.initAndValidate
        partialHelper.flipOffset(nodeNr);
        beagle.setPartials(partialHelper.getOffsetIndex(nodeNr), partials);
    }

    @Override
    protected void branchMatrix(int nr, int cat, double[] out) {
        System.arraycopy(javaMatrices[matrixHelper.getOffsetIndex(nr)], cat * S * S, out, 0, S * S);
    }

    @Override
    protected void nodePartials(int nr, double[] out) {
        beagle.getPartials(partialHelper.getOffsetIndex(nr), Beagle.NONE, out);
    }

    // ------------------------------------------------------------------ the tree likelihood on BEAGLE

    @Override
    protected int traverse(Node root) {
        // transition matrices of the changed branches (as TreeLikelihood.traverse decides and computes
        // them) and the queue of partials operations, then one updatePartials call
        int[] ops = {0};
        int update = traverseBeagle(root, ops);
        int opCount = ops[0];
        if (opCount > 0) {
            beagle.updatePartials(operations, opCount, Beagle.NONE);
        }
        if (update != Tree.IS_CLEAN) {
            for (int i = 0; i < nInternal; i++) scaleIndices[i] = scaleHelper.getOffsetIndex(i);
            beagle.resetScaleFactors(cumulativeScaleIndex);
            beagle.accumulateScaleFactors(scaleIndices, nInternal, cumulativeScaleIndex);
            int rootBuffer = partialHelper.getOffsetIndex(root.getNr());
            beagle.getPartials(rootBuffer, Beagle.NONE, rootPartials);
            double[] proportions = m_siteModel.getCategoryProportions(root);
            // each site's accumulated log scale factor: BEAGLE's site log likelihood under the uniform
            // root distribution (frequency buffer 0) minus the log of the same sum over the rescaled
            // root partials (getLogScaleFactors is not relied on)
            beagle.setCategoryWeights(0, proportions);
            beagle.calculateRootLogLikelihoods(new int[]{rootBuffer}, new int[]{0}, new int[]{0},
                    new int[]{cumulativeScaleIndex}, 1, sumLogLikelihood);
            beagle.getSiteLogLikelihoods(logScales);
            for (int k = 0; k < nPatterns; k++) {
                double sum = 0.0;
                for (int c = 0; c < nCat; c++) {
                    double s = 0.0;
                    for (int g = 0, v = (c * nPatterns + k) * S; g < S; g++, v++) s += uniform[g] * rootPartials[v];
                    sum += proportions[c] * s;
                }
                logScales[k] -= Math.log(sum);
            }
            double[] f = rootFrequenciesInput.get() != null
                    ? rootFrequenciesInput.get().getFreqs() : substitutionModel.getFrequencies();
            for (int k = 0; k < nPatterns; k++) {
                double sum = 0.0;
                for (int c = 0; c < nCat; c++) {
                    double s = 0.0;
                    for (int g = 0, v = (c * nPatterns + k) * S; g < S; g++, v++) s += f[g] * rootPartials[v];
                    sum += proportions[c] * s;
                }
                patternLogLikelihoods[k] = Math.log(sum) + logScales[k];
            }
        }
        return update;
    }

    /** Post-order: set changed matrices, queue the partials operations of nodes that need them. */
    private int traverseBeagle(Node node, int[] opCount) {
        int update = node.isDirty() | hasDirt;
        final int nr = node.getNr();
        if (!node.isRoot()) {
            double branchRate = branchRateModel.getRateForBranch(node);
            double branchTime = node.getLength() * branchRate;
            if (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nr]) {
                m_branchLengths[nr] = branchTime;
                matrixHelper.flipOffset(nr);
                int idx = matrixHelper.getOffsetIndex(nr);
                double[] all = javaMatrices[idx];
                if (all == null) all = javaMatrices[idx] = new double[nCat * S * S];
                for (int c = 0; c < nCat; c++) {
                    double rate = m_siteModel.getRateForCategory(c, node) * branchRate;
                    substitutionModel.getTransitionProbabilities(node, node.getParent().getHeight(), node.getHeight(),
                            rate, matrixScratch);
                    System.arraycopy(matrixScratch, 0, all, c * S * S, S * S);
                }
                beagle.setTransitionMatrix(idx, all, 0.0);
                update |= Tree.IS_DIRTY;
            }
        }
        if (!node.isLeaf()) {
            int u1 = traverseBeagle(node.getLeft(), opCount);
            int u2 = traverseBeagle(node.getRight(), opCount);
            if (u1 != Tree.IS_CLEAN || u2 != Tree.IS_CLEAN) {
                update |= (u1 | u2);
                partialHelper.flipOffset(nr);
                int si = nr - nTips;
                scaleHelper.flipOffset(si);
                int x = opCount[0] * Beagle.OPERATION_TUPLE_SIZE;
                operations[x] = partialHelper.getOffsetIndex(nr);
                operations[x + 1] = scaleHelper.getOffsetIndex(si);   // write this node's scale factors
                operations[x + 2] = Beagle.NONE;
                operations[x + 3] = partialHelper.getOffsetIndex(node.getLeft().getNr());
                operations[x + 4] = matrixHelper.getOffsetIndex(node.getLeft().getNr());
                operations[x + 5] = partialHelper.getOffsetIndex(node.getRight().getNr());
                operations[x + 6] = matrixHelper.getOffsetIndex(node.getRight().getNr());
                opCount[0]++;
            }
        }
        return update;
    }

    @Override
    public void store() {
        partialHelper.storeState();
        matrixHelper.storeState();
        scaleHelper.storeState();
        super.store();
    }

    @Override
    public void restore() {
        partialHelper.restoreState();
        matrixHelper.restoreState();
        scaleHelper.restoreState();
        super.restore();
    }

    /** Index i and its mirror for store/restore (BEAST BeagleTreeLikelihood's BufferIndexHelper). */
    private static final class BufferIndexHelper {
        private final int count;
        private int[] offsets, storedOffsets;

        BufferIndexHelper(int maxIndexValue, int minIndexValue) {
            if (minIndexValue != 0) throw new IllegalArgumentException("only every index mirrored");
            count = maxIndexValue;
            offsets = new int[count];
            storedOffsets = new int[count];
        }

        int getBufferCount() {
            return 2 * count;
        }

        void flipOffset(int i) {
            offsets[i] = count - offsets[i];
        }

        int getOffsetIndex(int i) {
            return offsets[i] + i;
        }

        void storeState() {
            System.arraycopy(offsets, 0, storedOffsets, 0, count);
        }

        void restoreState() {
            int[] tmp = storedOffsets;
            storedOffsets = offsets;
            offsets = tmp;
        }
    }
}
