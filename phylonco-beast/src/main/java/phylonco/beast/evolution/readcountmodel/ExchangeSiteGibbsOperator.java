package phylonco.beast.evolution.readcountmodel;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import mutablealignment.MutableAlignment;
import phylonco.beast.evolution.datatype.ReadCount;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Compound operator: narrow exchange + per-site joint Gibbs resample of all
 * leaf genotypes at K randomly chosen sites.
 *
 * <p>For each chosen site s, samples a new leaf-state vector g'_s jointly from
 * the conditional P(g_s | T', reads_s) using Felsenstein's pruning algorithm
 * with read-count likelihoods at the leaves. Because the joint per-site
 * proposal matches the per-site joint posterior conditional, the Hastings
 * ratio collapses cleanly and the effective acceptance probability for the
 * tree move is:</p>
 *
 * <pre>
 *   alpha = P(T')/P(T) * (validGP/validGPafter)
 *           * prod_{s in S} [P(reads_s | T') / P(reads_s | T)]
 * </pre>
 *
 * <p>i.e. the marginal likelihood ratio over resampled sites, with the
 * genotypes integrated out. Sites not in S still appear in BEAST's standard
 * likelihood ratio with their old genotypes, so K should be chosen high
 * enough that the tree move is judged primarily by integrated-out sites.</p>
 *
 * <p>The per-site joint Gibbs proposal density at a given configuration g
 * is computed as:</p>
 *
 * <pre>
 *   Q(g | T, reads_s) = P(g | T) * P(reads_s | g) / P(reads_s | T)
 * </pre>
 *
 * <p>where P(g | T) and P(reads_s | T) are computed by Felsenstein pruning
 * (leaves clamped to g, and leaves weighted by read-count likelihoods,
 * respectively).</p>
 */
@Description("Narrow exchange + per-site joint Gibbs resample at K sites. " +
        "Integrates leaf genotypes out of the tree-move acceptance ratio at the chosen sites.")
public class ExchangeSiteGibbsOperator extends TreeOperator {

    public Input<MutableAlignment> mutableAlignmentInput = new Input<>(
            "mutableAlignment",
            "mutable alignment containing genotype states",
            Input.Validate.REQUIRED);

    public Input<SiteModel.Base> siteModelInput = new Input<>(
            "siteModel",
            "site model containing substitution model",
            Input.Validate.REQUIRED);

    public Input<BranchRateModel.Base> branchRateModelInput = new Input<>(
            "branchRateModel",
            "branch rate model for relaxed clock (optional)",
            Input.Validate.OPTIONAL);

    public Input<LikelihoodReadCountModel> readCountModelInput = new Input<>(
            "readCountModel",
            "read count likelihood model",
            Input.Validate.REQUIRED);

    public Input<ReadCount> readCountInput = new Input<>(
            "readCount",
            "read count data",
            Input.Validate.REQUIRED);

    public Input<Integer> numSitesToResampleInput = new Input<>(
            "numSitesToResample",
            "number of sites to jointly Gibbs-resample per proposal (K). " +
                    "Use -1 (default) to resample all eligible sites — this collapses " +
                    "the Hastings ratio to the marginal-likelihood ratio with genotypes " +
                    "integrated out, which is the mathematically ideal tree-move " +
                    "acceptance criterion. For very large alignments, set a smaller K " +
                    "to trade mixing efficiency for per-call speed.",
            -1);

    public Input<Boolean> useEntropyWeightingInput = new Input<>(
            "useEntropyWeighting",
            "weight site-selection probability by per-site read-count entropy " +
                    "(sum over leaves of Shannon entropy of P(g | reads)); detailed balance " +
                    "preserved because weights depend only on read-count data. Default true.",
            true);

    private MutableAlignment alignment;
    private Tree tree;
    private SiteModel.Base siteModel;
    private SubstitutionModel substitutionModel;
    private BranchRateModel.Base branchRateModel;
    private LikelihoodReadCountModel readCountModel;
    private ReadCount readCount;

    private int numStates;
    private int numSites;
    private int numTaxa;
    private int numNodes;
    private int numCategories;
    private int K;
    private boolean useEntropyWeighting;

    // Per-site entropy weight: sum over leaves of Shannon entropy of P(g | reads_ls).
    // Function of read-count data only (computed once at init); preserves detailed
    // balance because forward and reverse site-selection probabilities agree.
    private double[] siteEntropy;

    // Pruning workspace.
    // partials[cat][node][state] is reused across calls, refilled per site.
    private double[][][] partials;
    private int[] sampledStates;         // for FFBS sampling
    private double[] transitionMatrix;
    private double[] logProbs;
    private double[] categoryLogProbs;
    private double[] rootFrequencies;

    // Coverages and index mappings (precomputed in initAndValidate).
    private int[][] coverages;            // [alignTaxonIdx][siteIdx]
    private int[] alignToRCIndex;         // alignment taxon idx -> ReadCount taxon idx
    private int[] nodeNrToTaxonIndex;     // tree node nr -> alignment taxon idx

    @Override
    public void initAndValidate() {
        alignment = mutableAlignmentInput.get();
        tree = treeInput.get();
        siteModel = siteModelInput.get();
        substitutionModel = siteModel.getSubstitutionModel();
        branchRateModel = branchRateModelInput.get();
        readCountModel = readCountModelInput.get();
        readCount = readCountInput.get();
        K = numSitesToResampleInput.get();
        useEntropyWeighting = useEntropyWeightingInput.get();

        if (K < 0) K = alignment.getSiteCount();  // -1 sentinel = all sites
        if (K < 1) {
            throw new IllegalArgumentException("numSitesToResample must be >= 1 or -1 (= all sites)");
        }

        numStates = alignment.getDataType().getStateCount();
        numSites = alignment.getSiteCount();
        numTaxa = alignment.getTaxonCount();
        numNodes = tree.getNodeCount();
        numCategories = siteModel.getCategoryCount();


        partials = new double[numCategories][numNodes][numStates];
        sampledStates = new int[numNodes];
        transitionMatrix = new double[numStates * numStates];
        logProbs = new double[numStates];
        categoryLogProbs = new double[numCategories];

        // Map alignment taxon names to ReadCount taxon indices.
        String[] rcTaxaNames = readCount.getTaxaNames();
        alignToRCIndex = new int[numTaxa];
        for (int i = 0; i < numTaxa; i++) {
            String alignTaxonName = alignment.getTaxaNames().get(i);
            alignToRCIndex[i] = i;
            for (int j = 0; j < rcTaxaNames.length; j++) {
                if (alignTaxonName.equals(rcTaxaNames[j].trim())) {
                    alignToRCIndex[i] = j;
                    break;
                }
            }
        }

        // Precompute coverages.
        coverages = new int[numTaxa][numSites];
        for (int i = 0; i < numTaxa; i++) {
            for (int j = 0; j < numSites; j++) {
                int[] counts = readCount.getReadCounts(alignToRCIndex[i], j);
                for (int k = 0; k < counts.length; k++) {
                    coverages[i][j] += counts[k];
                }
            }
        }

        // Tree node nr -> alignment taxon index.
        nodeNrToTaxonIndex = new int[numNodes];
        Arrays.fill(nodeNrToTaxonIndex, -1);
        for (Node leaf : tree.getExternalNodes()) {
            String taxonName = leaf.getID();
            int taxonIndex = alignment.getTaxonIndex(taxonName);
            nodeNrToTaxonIndex[leaf.getNr()] = taxonIndex;
        }

        // Per-site weight: Shannon entropy of the pooled nucleotide-count
        // distribution across all leaves. Captures both per-leaf uncertainty
        // and inter-leaf disagreement (e.g. subclonal SNVs where each cell's
        // reads are unambiguous but cells disagree on the dominant allele).
        // Pure read-count arithmetic; no model dependence, no MCMC-state
        // dependence; detailed balance preserved.
        //
        // A floor (a small fraction of max entropy log(numNucleotides)) is
        // added so monomorphic sites still get some pick weight — without it
        // they would never be Gibbs-resampled and could stay at their initial
        // genotypes forever if no other alignment operator runs.
        siteEntropy = new double[numSites];
        double entropyFloor = 0.0;
        for (int s = 0; s < numSites; s++) {
            int[] pooled = null;
            int total = 0;
            for (int t = 0; t < numTaxa; t++) {
                int[] counts = readCount.getReadCounts(alignToRCIndex[t], s);
                if (pooled == null) {
                    pooled = new int[counts.length];
                    // first site we see: set the floor to a small fraction
                    // of max possible entropy (log of alphabet size).
                    entropyFloor = 0.05 * Math.log(counts.length);
                }
                for (int n = 0; n < counts.length; n++) {
                    pooled[n] += counts[n];
                    total += counts[n];
                }
            }
            double h = 0.0;
            if (total > 0 && pooled != null) {
                for (int n = 0; n < pooled.length; n++) {
                    if (pooled[n] > 0) {
                        double p = (double) pooled[n] / total;
                        h -= p * Math.log(p);
                    }
                }
            }
            siteEntropy[s] = h + entropyFloor;
        }
    }

    @Override
    public double proposal() {
        final Tree tree = (Tree) InputUtil.get(treeInput, this);
        final int internalNodes = tree.getInternalNodeCount();

        if (internalNodes <= 1) {
            return Double.NEGATIVE_INFINITY;
        }

        // Pick grandparent with at least one non-leaf child.
        Node grandParent = tree.getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
        while (grandParent.getLeft().isLeaf() && grandParent.getRight().isLeaf()) {
            grandParent = tree.getNode(internalNodes + 1 + Randomizer.nextInt(internalNodes));
        }

        Node parentIndex = grandParent.getLeft();
        Node uncle = grandParent.getRight();
        if (parentIndex.getHeight() < uncle.getHeight()) {
            parentIndex = grandParent.getRight();
            uncle = grandParent.getLeft();
        }

        if (parentIndex.isLeaf()) {
            // dated-tip case
            return Double.NEGATIVE_INFINITY;
        }

        int validGP = 0;
        for (int n = internalNodes + 1; n < 1 + 2 * internalNodes; n++) {
            validGP += isg(tree.getNode(n));
        }
        final int c2 = sisg(parentIndex) + sisg(uncle);

        final Node i = Randomizer.nextBoolean() ? parentIndex.getLeft() : parentIndex.getRight();

        // Pick K eligible sites (weight > 0) without replacement.
        int[] weights = alignment.getWeights();
        int[] chosenSites = pickSites(weights, K);
        if (chosenSites == null) {
            return Double.NEGATIVE_INFINITY;
        }

        // Snapshot old genotypes and compute Q_R per site under current tree T.
        int[][] oldGenos = new int[chosenSites.length][];
        double logQR = 0.0;
        for (int s = 0; s < chosenSites.length; s++) {
            oldGenos[s] = alignment.getSiteValuesBySite(chosenSites[s]).clone();
            logQR += computeSiteLogDensity(chosenSites[s], oldGenos[s]);
        }

        // Narrow exchange T -> T'.
        exchangeNodes(i, uncle, parentIndex, grandParent);
        final int validGPafter = validGP - c2 + sisg(parentIndex) + sisg(uncle);

        // Sample new genotypes per site, compute Q_F under new tree T'.
        double logQF = 0.0;
        for (int s = 0; s < chosenSites.length; s++) {
            int[] newGenos = sampleSite(chosenSites[s]);
            logQF += computeSiteLogDensity(chosenSites[s], newGenos);
            alignment.setSiteValuesBySite(chosenSites[s], newGenos);
        }

        return Math.log((double) validGP / validGPafter) + logQR - logQF;
    }

    /**
     * Pick K distinct sites without replacement. Returns null if fewer than K
     * eligible sites exist.
     *
     * When entropy weighting is enabled, site s is sampled with probability
     * proportional to patternWeight[s] * siteEntropy[s]. Both factors depend
     * only on data, so forward and reverse selection probabilities for any
     * given set agree and cancel in the Hastings ratio.
     */
    private int[] pickSites(int[] weights, int K) {
        if (!useEntropyWeighting) {
            // Uniform Fisher-Yates over patternWeight > 0 sites.
            int eligibleCount = 0;
            for (int w : weights) if (w > 0) eligibleCount++;
            if (eligibleCount < K) return null;
            int[] eligible = new int[eligibleCount];
            int idx = 0;
            for (int s = 0; s < weights.length; s++) {
                if (weights[s] > 0) eligible[idx++] = s;
            }
            int[] picked = new int[K];
            for (int k = 0; k < K; k++) {
                int swap = k + Randomizer.nextInt(eligibleCount - k);
                picked[k] = eligible[swap];
                eligible[swap] = eligible[k];
            }
            return picked;
        }

        // Entropy-weighted sampling without replacement.
        double[] effective = new double[numSites];
        int eligibleCount = 0;
        double totalWeight = 0.0;
        for (int s = 0; s < numSites; s++) {
            if (weights[s] <= 0) continue;
            double w = weights[s] * siteEntropy[s];
            if (w > 0.0) {
                effective[s] = w;
                totalWeight += w;
                eligibleCount++;
            }
        }
        if (eligibleCount < K || totalWeight <= 0.0) return null;

        int[] picked = new int[K];
        for (int k = 0; k < K; k++) {
            double r = Randomizer.nextDouble() * totalWeight;
            double cum = 0.0;
            int chosen = -1;
            for (int s = 0; s < numSites; s++) {
                if (effective[s] > 0.0) {
                    cum += effective[s];
                    if (r < cum) { chosen = s; break; }
                }
            }
            if (chosen < 0) {
                for (int s = numSites - 1; s >= 0; s--) {
                    if (effective[s] > 0.0) { chosen = s; break; }
                }
            }
            picked[k] = chosen;
            totalWeight -= effective[chosen];
            effective[chosen] = 0.0;
        }
        return picked;
    }

    /**
     * Returns log P(g | T, reads_s) where g is the given leaf genotype vector
     * (indexed by alignment taxon index) at site s, T is the tree's current
     * topology/branch-lengths, and reads_s are the read counts at site s.
     *
     * Computed as log P(g | T) + log P(reads_s | g) - log P(reads_s | T)
     * via two Felsenstein pruning passes (leaves clamped to g; leaves
     * weighted by read-count likelihoods).
     */
    private double computeSiteLogDensity(int siteIdx, int[] genotypes) {
        // Pass 1: leaves clamped to g -> log P(g | T).
        for (Node leaf : tree.getExternalNodes()) {
            setLeafPartialClamped(leaf, genotypes);
        }
        for (int cat = 0; cat < numCategories; cat++) {
            computeInternalPartials(tree.getRoot(), cat);
        }
        double logP_g_T = computeRootLogMarginal();

        // P(reads_s | g) = product over leaves of P(reads_l | g_l) (read-count model).
        double logP_reads_g = 0.0;
        for (Node leaf : tree.getExternalNodes()) {
            int taxonIdx = nodeNrToTaxonIndex[leaf.getNr()];
            int rcIdx = alignToRCIndex[taxonIdx];
            int[] counts = readCount.getReadCounts(rcIdx, siteIdx);
            int coverage = coverages[taxonIdx][siteIdx];
            logP_reads_g += readCountModel.logLiklihoodRC(genotypes[taxonIdx], counts, coverage, rcIdx);
        }

        // Pass 2: leaves weighted by read-count likelihoods -> log P(reads_s | T).
        for (Node leaf : tree.getExternalNodes()) {
            setLeafPartialReadCount(leaf, siteIdx);
        }
        for (int cat = 0; cat < numCategories; cat++) {
            computeInternalPartials(tree.getRoot(), cat);
        }
        double logP_reads_T = computeRootLogMarginal();

        return logP_g_T + logP_reads_g - logP_reads_T;
    }

    /**
     * Sample joint leaf genotypes at one site from P(g | T, reads_s) using
     * forward-filter / backward-sample. Returns the sampled g vector
     * (indexed by alignment taxon index).
     */
    private int[] sampleSite(int siteIdx) {
        for (Node leaf : tree.getExternalNodes()) {
            setLeafPartialReadCount(leaf, siteIdx);
        }
        for (int cat = 0; cat < numCategories; cat++) {
            computeInternalPartials(tree.getRoot(), cat);
        }

        int category = sampleCategoryAndRootState();
        sampleDescendantStates(tree.getRoot(), category);

        int[] genos = new int[numTaxa];
        for (Node leaf : tree.getExternalNodes()) {
            int taxonIdx = nodeNrToTaxonIndex[leaf.getNr()];
            genos[taxonIdx] = sampledStates[leaf.getNr()];
        }
        return genos;
    }

    private void setLeafPartialClamped(Node leaf, int[] genotypes) {
        int nodeNr = leaf.getNr();
        int taxonIdx = nodeNrToTaxonIndex[nodeNr];
        int clampedState = genotypes[taxonIdx];
        for (int cat = 0; cat < numCategories; cat++) {
            for (int g = 0; g < numStates; g++) {
                partials[cat][nodeNr][g] = (g == clampedState) ? 1.0 : 0.0;
            }
        }
    }

    private void setLeafPartialReadCount(Node leaf, int siteIdx) {
        int nodeNr = leaf.getNr();
        int taxonIdx = nodeNrToTaxonIndex[nodeNr];
        int rcIdx = alignToRCIndex[taxonIdx];
        int[] counts = readCount.getReadCounts(rcIdx, siteIdx);
        int coverage = coverages[taxonIdx][siteIdx];
        for (int g = 0; g < numStates; g++) {
            double logP = readCountModel.logLiklihoodRC(g, counts, coverage, rcIdx);
            double p = Math.exp(logP);
            for (int cat = 0; cat < numCategories; cat++) {
                partials[cat][nodeNr][g] = p;
            }
        }
    }

    private void computeInternalPartials(Node node, int category) {
        if (node.isLeaf()) return;

        for (int c = 0; c < node.getChildCount(); c++) {
            computeInternalPartials(node.getChild(c), category);
        }

        int nodeNr = node.getNr();
        Arrays.fill(partials[category][nodeNr], 1.0);

        for (int c = 0; c < node.getChildCount(); c++) {
            Node child = node.getChild(c);
            int childNr = child.getNr();
            getTransitionMatrix(child, category, transitionMatrix);

            for (int parentState = 0; parentState < numStates; parentState++) {
                double sum = 0.0;
                for (int childState = 0; childState < numStates; childState++) {
                    double trans = transitionMatrix[parentState * numStates + childState];
                    sum += trans * partials[category][childNr][childState];
                }
                partials[category][nodeNr][parentState] *= sum;
            }
        }
    }

    private void getTransitionMatrix(Node node, int category, double[] matrix) {
        double branchRate = (branchRateModel != null) ? branchRateModel.getRateForBranch(node) : 1.0;
        double rate = branchRate * siteModel.getRateForCategory(category, node);
        substitutionModel.getTransitionProbabilities(
                node,
                node.getParent().getHeight(),
                node.getHeight(),
                rate,
                matrix);
    }

    /**
     * Returns log [ sum_cat prop[cat] * sum_g pi[g] * partials[cat][root][g] ]
     * given partials are already built.
     */
    private double computeRootLogMarginal() {
        Node root = tree.getRoot();
        int rootNr = root.getNr();
        rootFrequencies = substitutionModel.getFrequencies();
        double[] proportions = siteModel.getCategoryProportions(root);

        double max = Double.NEGATIVE_INFINITY;
        double[] catLogProbs = new double[numCategories];
        for (int cat = 0; cat < numCategories; cat++) {
            double marginal = 0.0;
            for (int g = 0; g < numStates; g++) {
                marginal += rootFrequencies[g] * partials[cat][rootNr][g];
            }
            catLogProbs[cat] = Math.log(proportions[cat]) + Math.log(marginal);
            if (catLogProbs[cat] > max) max = catLogProbs[cat];
        }
        if (max == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;
        double sum = 0.0;
        for (int cat = 0; cat < numCategories; cat++) {
            sum += Math.exp(catLogProbs[cat] - max);
        }
        return max + Math.log(sum);
    }

    /**
     * Sample (category, root state) jointly and write sampledStates[root] = root state.
     * Returns the sampled category.
     */
    private int sampleCategoryAndRootState() {
        Node root = tree.getRoot();
        int rootNr = root.getNr();
        rootFrequencies = substitutionModel.getFrequencies();
        double[] proportions = siteModel.getCategoryProportions(root);

        double maxCatLogProb = Double.NEGATIVE_INFINITY;
        for (int cat = 0; cat < numCategories; cat++) {
            double marginal = 0.0;
            for (int g = 0; g < numStates; g++) {
                marginal += rootFrequencies[g] * partials[cat][rootNr][g];
            }
            categoryLogProbs[cat] = Math.log(proportions[cat]) + Math.log(marginal);
            if (categoryLogProbs[cat] > maxCatLogProb) maxCatLogProb = categoryLogProbs[cat];
        }
        int sampledCategory = sampleFromLogProbs(categoryLogProbs, numCategories, maxCatLogProb);

        double maxLogProb = Double.NEGATIVE_INFINITY;
        for (int g = 0; g < numStates; g++) {
            logProbs[g] = Math.log(rootFrequencies[g]) + Math.log(partials[sampledCategory][rootNr][g]);
            if (logProbs[g] > maxLogProb) maxLogProb = logProbs[g];
        }
        sampledStates[rootNr] = sampleFromLogProbs(logProbs, numStates, maxLogProb);
        return sampledCategory;
    }

    /**
     * Pre-order descendant sampling: each child state sampled from
     * P(child=g | parent state, partial[child][g]).
     */
    private void sampleDescendantStates(Node node, int category) {
        for (int c = 0; c < node.getChildCount(); c++) {
            Node child = node.getChild(c);
            int childNr = child.getNr();
            int parentState = sampledStates[node.getNr()];

            getTransitionMatrix(child, category, transitionMatrix);

            double maxLogProb = Double.NEGATIVE_INFINITY;
            for (int g = 0; g < numStates; g++) {
                double trans = transitionMatrix[parentState * numStates + g];
                logProbs[g] = Math.log(trans) + Math.log(partials[category][childNr][g]);
                if (logProbs[g] > maxLogProb) maxLogProb = logProbs[g];
            }
            sampledStates[childNr] = sampleFromLogProbs(logProbs, numStates, maxLogProb);

            if (!child.isLeaf()) {
                sampleDescendantStates(child, category);
            }
        }
    }

    private int sampleFromLogProbs(double[] logP, int size, double maxLog) {
        double sumExp = 0.0;
        double[] probs = new double[size];
        for (int i = 0; i < size; i++) {
            probs[i] = Math.exp(logP[i] - maxLog);
            sumExp += probs[i];
        }
        double rand = Randomizer.nextDouble() * sumExp;
        double cum = 0.0;
        for (int i = 0; i < size; i++) {
            cum += probs[i];
            if (rand < cum) return i;
        }
        return size - 1;
    }

    private void exchangeNodes(Node i, Node j, Node p, Node jP) {
        replace(p, i, j);
        replace(jP, j, i);
    }

    private int isg(Node n) {
        return (n.getLeft().isLeaf() && n.getRight().isLeaf()) ? 0 : 1;
    }

    private int sisg(Node n) {
        return n.isLeaf() ? 0 : isg(n);
    }
}
